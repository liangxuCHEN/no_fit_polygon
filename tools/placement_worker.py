# -*- coding: utf-8 -*-
import json
from nfp_utls import almost_equal, rotate_polygon, get_polygon_bounds, polygon_area
import copy
import pyclipper


class PlacementWorker():
    def __init__(self, bin_polygon, paths, ids, rotations, config, nfp_cache):
        self.bin_polygon = bin_polygon
        self.paths = copy.deepcopy(paths)
        self.ids = ids       # 图形原来的ID顺序
        self.rotations = rotations
        self.config = config
        self.nfpCache = nfp_cache or {}

    def place_paths(self):
        # 排列图形
        if self.bin_polygon is None:
            return None

        # rotate paths by given rotation
        rotated = list()
        for i in range(0, len(self.paths)):
            r = rotate_polygon(self.paths[i][1]['points'], self.paths[i][2])
            r['rotation'] = self.paths[i][2]
            r['source'] = self.paths[i][1]['p_id']
            r['p_id'] = self.paths[i][0]
            rotated.append(r)

        paths = rotated
        # 保存所有转移数据
        all_placements = list()
        # 基因组的适应值
        fitness = 0
        bin_area = abs(polygon_area(self.bin_polygon['points']))
        min_width = None
        while len(paths) > 0:
            placed = list()
            placements = list()
            # add 1 for each new bin opened (lower fitness is better)
            fitness += 1
            for i in range(0, len(paths)):
                path = paths[i]
                # 图形的坐标
                key = json.dumps({
                    'A': '-1',
                    'B': path['p_id'],
                    'inside': True,
                    'A_rotation': 0,
                    'B_rotation': path['rotation']
                })

                binNfp = self.nfpCache.get(key)
                if binNfp is None or len(binNfp) == 0:
                    continue

                # part unplaceable, skip
                error = False

                # ensure all necessary NFPs exist
                for p in placed:
                    key = json.dumps({
                        'A': p['p_id'],
                        'B': path['p_id'],
                        'inside': False,
                        'A_rotation': p['rotation'],
                        'B_rotation': path['rotation']
                    })
                    nfp = self.nfpCache.get(key)
                    if nfp is None:
                        error = True
                        break

                # part unplaceable, skip
                if error:
                    continue

                position = None
                if len(placed) == 0:
                    for j in range(0, len(binNfp)):
                        for k in range(0, len(binNfp[j])):
                            if position is None or (binNfp[j][k]['x']-path['points'][0]['x'] < position['x']):
                                position = {
                                    'x': binNfp[j][k]['x'] - path['points'][0]['x'],
                                    'y': binNfp[j][k]['y'] - path['points'][0]['y'],
                                    'p_id': path['p_id'],
                                    'rotation': path['rotation']
                                }

                    placements.append(position)
                    placed.append(path)
                    continue

                clipper_bin_nfp = list()
                for j in range(0, len(binNfp)):
                    clipper_bin_nfp.append([[p['x'], p['y']] for p in binNfp[j]])

                clipper = pyclipper.Pyclipper()

                for j in range(0, len(placed)):
                    p = placed[j]
                    key = json.dumps({
                        'A': p['p_id'],
                        'B': path['p_id'],
                        'inside': False,
                        'A_rotation': p['rotation'],
                        'B_rotation': path['rotation']
                    })
                    nfp = self.nfpCache.get(key)

                    if nfp is None:
                        continue
                    for k in range(0, len(nfp)):
                        clone = [[np['x'] + placements[j]['x'], np['y'] + placements[j]['y']] for np in nfp[k]]
                        clone = pyclipper.CleanPolygon(clone)
                        if len(clone) > 2:
                            clipper.AddPath(clone, pyclipper.PT_SUBJECT, True)
                combine_nfp = clipper.Execute(pyclipper.CT_UNION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
                if len(combine_nfp) == 0:
                    continue

                clipper = pyclipper.Pyclipper()
                clipper.AddPaths(combine_nfp, pyclipper.PT_CLIP, True)
                try:
                    clipper.AddPaths(clipper_bin_nfp, pyclipper.PT_SUBJECT, True)
                except:
                    print u'图形坐标出错', clipper_bin_nfp

                # choose placement that results in the smallest bounding box
                finalNfp = clipper.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
                if len(finalNfp) == 0:
                    continue
                finalNfp = pyclipper.CleanPolygons(finalNfp)

                for j in range(len(finalNfp)-1, -1, -1):
                    if len(finalNfp[j]) < 3:
                        finalNfp.pop(j)
                if len(finalNfp) == 0:
                    continue

                finalNfp = [[{'x': p[0], 'y': p[1]}for p in polygon] for polygon in finalNfp]

                min_width = None
                min_area = None
                min_x = None

                for nf in finalNfp:

                    if abs(polygon_area(nf)) < 2:
                        continue

                    for p_nf in nf:
                        # 生成nfp多边形
                        all_points = list()
                        for m in range(0, len(placed)):
                            for p in placed[m]['points']:
                                all_points.append({
                                    'x': p['x']+placements[m]['x'],
                                    'y': p['y']+placements[m]['y']
                                })
                        # path 坐标
                        shift_vector = {
                            'x': p_nf['x'] - path['points'][0]['x'],
                            'y': p_nf['y'] - path['points'][0]['y'],
                            'p_id': path['p_id'],
                            'rotation': path['rotation'],
                        }

                        # 找新坐标后的最小矩形
                        for m in range(0, len(path['points'])):
                            all_points.append({
                                'x': path['points'][m]['x'] + shift_vector['x'],
                                'y': path['points'][m]['y'] + shift_vector['y']
                            })

                        rect_bounds = get_polygon_bounds(all_points)
                        # weigh width more, to help compress in direction of gravity
                        area = rect_bounds['width'] * 2 + rect_bounds['height']

                        if (min_area is None or area < min_area or almost_equal(min_area, area)) and (
                                        min_x is None or shift_vector['x'] <= min_x):
                            min_area = area
                            min_width = rect_bounds['width']
                            position = shift_vector
                            min_x = shift_vector['x']

                if position:
                    placed.append(path)
                    placements.append(position)

            if min_width:
                fitness += min_width / bin_area

            for p in placed:
                p_id = paths.index(p)
                if p_id >= 0:
                    paths.pop(p_id)

            if placements and len(placements) > 0:
                all_placements.append(placements)

            else:
                # something went wrong
                break

        fitness += 2 * len(paths)

        return {'placements': all_placements, 'fitness': fitness, 'paths': paths, 'area': bin_area}
