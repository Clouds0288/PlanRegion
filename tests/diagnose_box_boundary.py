"""An exact-oracle counterexample: radial box splitting can stall near a face."""
import json
from pathlib import Path

import numpy as np

from tests.benchmark_boundary_search import clip_boxes, coverage, inner_scales, write_json
from region import GEOMETRY_TOL


def main():
    # This is the entire true domain, not an approximation or a solver result.
    inner_points = np.array([[1., .3, 1.], [.3, 1., 1.]])
    outer_vertices = np.ones((1, 3))
    history = []
    for index in range(100):
        gap, residual = coverage(outer_vertices, inner_points)
        if np.all(residual <= GEOMETRY_TOL):
            break
        point = outer_vertices[int(np.argmax(gap))]
        # Exact ray maximum of the known two-box domain.
        theta_upper = inner_scales(point[None, :], inner_points)[0]
        before = outer_vertices.copy()
        outer_vertices = clip_boxes(outer_vertices, point, theta_upper, bounds=np.ones(3))
        history.append(dict(iteration=index+1, direction=point, theta_upper=theta_upper,
                            max_coverage_gap=float(gap.max()), unchanged=np.array_equal(before, outer_vertices)))
        if history[-1]['unchanged']:
            break
    gap, residual = coverage(outer_vertices, inner_points)
    # One valid zero-coordinate query removes the spurious third OR branch.
    repaired = clip_boxes(outer_vertices, np.array([1., 1., 0.]), .3, bounds=np.ones(3))
    repaired_gap, repaired_residual = coverage(repaired, inner_points)
    result = dict(domain='[0,(1,.3,1)] union [0,(.3,1,1)]', inner_is_exact=True,
                  oracle_is_exact=True, history=history, outer_vertices=outer_vertices,
                  max_coverage_gap=float(gap.max()), max_covering_residual=float(residual.max()),
                  certified=bool(np.all(residual <= GEOMETRY_TOL)),
                  after_plane_query=dict(outer_vertices=repaired, max_coverage_gap=float(repaired_gap.max()),
                      certified=bool(np.all(repaired_residual <= GEOMETRY_TOL))))
    write_json(Path('results/case33_boundary_search/20260924/synthetic_face_diagnostic.json'), result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('history', 'outer_vertices', 'after_plane_query')}))
    print('iterations',len(history),'unchanged',history[-1]['unchanged'], 'plane_query_certified',result['after_plane_query']['certified'])
    assert not result['certified'] and result['after_plane_query']['certified']


if __name__ == '__main__':
    main()
