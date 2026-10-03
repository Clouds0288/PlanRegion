"""Check docs/notation.md code bindings without importing models or solvers."""
import ast
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
CODE_LINK = re.compile(r"\[([^\]\n]+)\]\(\.\./([^\s)]+\.py)\)")
QUALIFIED_NAME = re.compile(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*(?::[A-Za-z_]\w*)?")
DETAILS = {'param', 'local'}   # 形参与函数内局部量是实现细节，不登记


class DefinitionIndex(ast.NodeVisitor):
    """Collect definitions, self fields and explicitly authored record keys, with the kind of each name:
    module, definition, attribute, param, local or key. Every name bound inside a function, nested
    functions and classes included, is local."""

    def __init__(self, source):
        self.kinds = {}
        self.scope, self.frames = [], []
        self.class_scope = None
        self.visit(ast.parse(source))

    @property
    def names(self):
        return set(self.kinds)

    def record(self, name, kind):
        self.kinds.setdefault(name, set()).add(kind)

    def add(self, name, kind):
        self.record('.'.join([*self.scope, name]), kind)

    def define(self, name):
        self.add(name, 'local' if self.frames and self.frames[-1] == 'function' else 'definition')

    def visit_ClassDef(self, node):
        self.define(node.name)
        previous = self.class_scope
        self.scope.append(node.name)
        self.frames.append('class')
        self.class_scope = '.'.join(self.scope)
        for statement in node.body:
            self.visit(statement)
        self.scope.pop()
        self.frames.pop()
        self.class_scope = previous

    def visit_FunctionDef(self, node):
        self.define(node.name)
        self.scope.append(node.name)
        self.frames.append('function')
        args = node.args
        for arg in [*args.posonlyargs, *args.args, *args.kwonlyargs, args.vararg, args.kwarg]:
            if arg is not None:
                self.add(arg.arg, 'param')
        for statement in node.body:
            self.visit(statement)
        self.scope.pop()
        self.frames.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Store):
            self.add(node.id, 'module' if not self.frames else 'attribute' if self.frames[-1] == 'class' else 'local')

    def visit_Attribute(self, node):
        if (self.class_scope and isinstance(node.ctx, ast.Store)
                and isinstance(node.value, ast.Name) and node.value.id == 'self'):
            self.record(f'{self.class_scope}.{node.attr}', 'attribute')
        self.generic_visit(node)

    def dict_key(self, key):
        if self.scope:
            self.record(f'{".".join(self.scope)}:{key}', 'key')

    def visit_Dict(self, node):
        for key in node.keys:
            if isinstance(key, ast.Constant) and isinstance(key.value, str):
                self.dict_key(key.value)
        self.generic_visit(node)

    def visit_Call(self, node):
        """dict(...)、SimpleNamespace(...) 与 .update(...) 的关键字是记录的键。"""
        is_record = isinstance(node.func, ast.Name) and node.func.id in ('dict', 'SimpleNamespace')
        is_update = isinstance(node.func, ast.Attribute) and node.func.attr == 'update'
        if is_record or is_update:
            for keyword in node.keywords:
                if keyword.arg is not None:
                    self.dict_key(keyword.arg)
        self.generic_visit(node)


class NotationContractTests(unittest.TestCase):
    def test_documented_names_exist_in_their_declared_scope(self):
        document = (ROOT / 'docs' / 'notation.md').read_text(encoding='utf-8')
        references = CODE_LINK.findall(document)
        self.assertTrue(references, 'No code definitions registered in docs/notation.md')
        indexes = {}
        for name, relative in references:
            with self.subTest(file=relative, name=name):
                self.assertIsNotNone(QUALIFIED_NAME.fullmatch(name),
                                     f'Invalid notation reference: {name}')
                path = (ROOT / relative).resolve()
                self.assertTrue(path.is_relative_to(ROOT), 'Reference leaves the repository')
                self.assertTrue(path.is_file(), f'Missing source file: {relative}')
                if relative not in indexes:
                    indexes[relative] = DefinitionIndex(path.read_text(encoding='utf-8-sig')).kinds
                kinds = indexes[relative].get(name)
                self.assertIsNotNone(kinds, f'Notation contract drift: {relative} no longer defines {name}. '
                                            'Preserve the registered name or document an explicit migration.')
                self.assertFalse(kinds <= DETAILS, f'{relative}: {name} is a function parameter or local; '
                                                   'implementation details are not registered.')
                self.assertFalse(re.split('[.:]', name)[-1].startswith('_'),
                                 f'{relative}: {name} is private; implementation details are not registered.')

    def test_reads_and_other_scopes_do_not_mask_a_rename(self):
        source = '''
class Network:
    reactance: object
    def loads(self, power):
        return dict(p=power)
class Unrelated:
    reactance: object
    def loads(self, power):
        return dict(p=power)
def consumer(network, answer):
    return network.reactance, answer['p']
'''
        original = DefinitionIndex(source).names
        protected = {'Network.reactance', 'Network.loads.power', 'Network.loads:p'}
        self.assertTrue(protected <= original)
        # Simulate drift in memory, leaving working files untouched.
        tree = ast.parse(source)
        network = tree.body[0]
        network.body[0].target.id = 'x'
        loads = network.body[1]
        loads.args.args[1].arg = 'demand'
        loads.body[0].value.keywords[0].arg = 'load'
        renamed = DefinitionIndex(ast.unparse(tree)).names
        self.assertTrue({'Network.x', 'Network.loads.demand', 'Network.loads:load'} <= renamed)
        self.assertTrue(protected.isdisjoint(renamed))
        self.assertIn('Unrelated.reactance', renamed)

    def test_kinds_separate_interface_from_implementation_details(self):
        source = '''
LIMIT = 1.
class Model:
    scale = 2.
    def solve(self, power):
        state = power*self.scale
        self.state = state
        def helper():
            return state
        return SimpleNamespace(state=helper()), dict(p=power)
'''
        kinds = DefinitionIndex(source).kinds
        self.assertEqual(kinds['LIMIT'], {'module'})
        self.assertEqual(kinds['Model.solve'], {'definition'})
        self.assertEqual(kinds['Model.scale'], {'attribute'})
        self.assertEqual(kinds['Model.state'], {'attribute'})
        self.assertEqual(kinds['Model.solve.power'], {'param'})
        self.assertEqual(kinds['Model.solve.state'], {'local'})
        self.assertEqual(kinds['Model.solve.helper'], {'local'})
        self.assertEqual(kinds['Model.solve:state'], {'key'})
        self.assertEqual(kinds['Model.solve:p'], {'key'})


if __name__ == '__main__':
    unittest.main()
