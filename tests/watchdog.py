"""仅供测试的外部进程看门狗；正式 Notebook 不导入、不报告检查点。"""
import multiprocessing as mp  # 独立工作进程可在硬超时时被终止。
from pathlib import Path  # 固定被测 Notebook 的绝对路径。
from time import perf_counter  # 测试时限包含进程启动和导入。
import traceback  # 将被测代码的真实异常传回测试进程。
import nbformat  # 从正式 Notebook 加载流程定义。
from threadpoolctl import threadpool_limits  # 保持测试矩阵运算单线程。


def _execute(notebook, function, kwargs, channel):  # spawn 入口必须位于可导入的测试模块。
    try:  # 进程边界统一传回结果或异常。
        namespace = {}  # 被测流程共享自己的独立命名空间。
        for cell in nbformat.read(notebook,as_version=4).cells:  # 只加载导入和函数，不启动整个实验。
            tags = cell.metadata.get('tags',[])  # 与测试函数加载器采用同一规则。
            if cell.cell_type=='code' and ('imports' in tags or any(t.endswith('-flow') for t in tags)):  # 跳过执行和绘图单元。
                exec(compile(cell.source,notebook,'exec'),namespace)  # 不向正式函数注入 deadline 或报告接口。
        with threadpool_limits(limits=1):  # 与主流程使用相同数值线程设置。
            value = namespace[function](**kwargs)  # 正式函数无需知道自己受测试看门狗保护。
        channel.send(dict(value=value,termination='returned',error=None))  # 只传回本次调用的最终结果。
    except Exception:  # 异常不能被测试工具伪装成物理不可行。
        channel.send(dict(value=None,termination='error',error=traceback.format_exc()))  # 保存完整堆栈。
    finally:  # 关闭工作进程自己的通信端。
        channel.close()  # 超时强制终止时由操作系统回收。


def run_guarded(function, kwargs, seconds, *, notebook='main.ipynb'):  # 从测试侧限制一次函数调用的总时间。
    context = mp.get_context('spawn')  # 不依赖 Notebook 内核中的可序列化函数对象。
    parent, child = context.Pipe(duplex=False)  # 只由被测工作进程传回结果。
    start = perf_counter()  # 计时覆盖创建、导入和执行。
    process = context.Process(target=_execute,args=(str(Path(notebook).resolve()),function,kwargs,child))  # 仅创建当前测试拥有的进程。
    process.start()  # 无须主流程提供任何看门狗钩子。
    child.close()  # 父进程关闭不用的发送端。
    result = dict(value=None,termination='hard_timeout',error=None)  # 超时是测试未完成，不是模型不可行。
    try:  # 正常返回、异常和用户中断均需清理工作进程。
        if parent.poll(max(0.,seconds-(perf_counter()-start))):  # 等待最终消息，启动时间也计入上限。
            try:  # 工作进程崩溃可能只关闭通道而无结果。
                result = parent.recv()  # 收到一次完整结果后立即停止等待。
            except EOFError:  # 区分异常退出与正常返回。
                result['termination'] = 'worker_exit'  # 不重试失败测试。
    finally:  # 只回收本次测试创建的进程。
        process.join(.2)  # 留出正常退出的短暂资源清理时间。
        if process.is_alive():  # 硬时限内未结束时才强制终止。
            process.terminate()  # 不扫描或终止其他求解任务。
        process.join()  # 回收当前进程句柄。
        parent.close()  # 关闭父进程通信端。
    return dict(**result,wall_seconds=perf_counter()-start)  # 测试自行解释结果，不写正式实验文件。
