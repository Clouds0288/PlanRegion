"""主线三维并集测度：只积分暴露面，保持重叠和内部空隙。"""
import numpy as np


def _clip_face(points, constant, coefficient):
    """有序面片裁剪：constant+coefficient@point>=0，不重新求凸包。"""
    values = (constant+points@coefficient)/np.linalg.norm(coefficient)
    values[np.abs(values) <= 1e-11] = 0.
    clipped = []
    for index in range(len(points)):
        previous = index-1
        if values[previous]*values[index] < 0.:
            clipped.append(points[previous]+(points[index]-points[previous])
                           *values[previous]/(values[previous]-values[index]))
        if values[index] >= 0.:
            clipped.append(points[index])
    return np.asarray(clipped).reshape(-1, 2)


def union_volume(polytopes):
    """散度定理：只积分并集暴露的边界面，重叠面和内部面不重复计入。"""
    from scipy.spatial import ConvexHull
    from shapely import set_precision
    from shapely.geometry import Polygon
    from region import polytope_volume, halfspaces
    polytopes = [np.asarray(p) for p in polytopes if polytope_volume(p) > 1e-15]
    equations = [halfspaces(p) for p in polytopes]
    volume = 0.
    for i, poly in enumerate(polytopes):
        coordinates = np.linalg.svd(poly-poly[0], full_matrices=False)[0]
        for simplex in ConvexHull(coordinates, qhull_options='Qx').simplices:
            triangle = poly[simplex]
            a, b = triangle[1:]-triangle[0]
            normal = np.cross(a, b)
            area = np.linalg.norm(normal)
            if area == 0.:
                continue
            normal /= area
            if normal@(triangle[0]-poly.mean(axis=0)) < 0.:
                normal = -normal
            offset = -normal@triangle[0]
            basis = np.array([a/np.linalg.norm(a), np.cross(normal, a/np.linalg.norm(a))])
            face = (triangle-triangle[0])@basis.T
            exposed = Polygon(face)
            for j, eq in enumerate(equations):
                if i == j:
                    continue
                # 共面且同向的重叠外表面由索引较小的网架计入一次。
                coincident = ((np.linalg.norm(eq[:, :3]-normal, axis=1) < 1e-10)
                              & (np.abs(eq[:, 3]-offset) < 1e-10))
                if j > i and coincident.any():
                    continue
                values = triangle@eq[:, :3].T+eq[:, 3]
                if np.any(values.min(axis=0) > 1e-10):
                    continue
                if np.all(values <= 1e-10):
                    exposed = Polygon()
                    break
                coefficients = eq[:, :3]@basis.T
                constants = eq[:, :3]@triangle[0]+eq[:, 3]
                parallel = np.linalg.norm(coefficients, axis=1) < 1e-10
                if np.any(constants[parallel] > 1e-10):
                    continue
                overlap = face
                for coefficient, constant in zip(coefficients[~parallel], constants[~parallel]):
                    overlap = _clip_face(overlap, -constant, -coefficient)
                    if not len(overlap):
                        break
                if len(overlap) >= 3:
                    # 共边点统一到面内裁剪精度，避免微小共线偏差使重叠面漏扣。
                    exposed = set_precision(exposed, 1e-11).difference(set_precision(Polygon(overlap), 1e-11))
                    if exposed.is_empty:
                        break
            volume -= offset*exposed.area/3.
    return volume
