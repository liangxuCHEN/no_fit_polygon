# -*- coding: utf-8 -*-
from datetime import datetime
from tools import placement_worker, nfp_utls, input_utls
import math
import json
import random
import copy
from Polygon import Polygon
from Polygon.IO import writeSVG
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pyclipper


TOLERANCE = 0.0001  # smaller than this, two points are considered equal
DISCRETIZE = 4  # the number of segments in which arcs must be subdivided
ROTATIONS = [0, 90, 180, 270]  # the possible rotations to try


class Nester:
    def __init__(self, container=None, shapes=None):

        """Nester([container,shapes]): Creates a nester object with a container
           shape and a list of other shapes to nest into it. Container and
           shapes must be Part.Faces.
           Typical workflow:
           n = Nester() # creates the nester
           n.addContainer(object) # adds a doc object as the container
           n.addObjects(objects) # adds a list of doc objects as shapes
           n.run() # runs the nesting
           n.show() # creates a preview (compound) of the results
           """

        self.container = container
        self.bin_bounds = None
        self.shapes = shapes
        self.results = []  # storage for the different results
        self.indexedFaces = None
        self.setCounter = None  # optionally define a setCounter(value) function where value is a %
        self.nfpCache = {}
        self.parts = None
        self.config = {
            'curveTolerance': 0.3,  # 允许的最大误差转换贝济耶和圆弧线段。在SVG的单位。更小的公差将需要更长的时间来计算
            'spacing': 3,
            'rotations': 4,
            'populationSize': 10,
            'mutationRate': 10,
            'useHoles': True,
            'exploreConcave': False  # 寻找凹面
        }

        self.working = False

        self.GA = None
        self.best = None
        self.worker = None


    def addObjects(self, objects):
        """addObjects(objects): adds polygon objects to the nester"""
        if not isinstance(objects, list):
            objects = [objects]
        if not self.shapes:
            self.shapes = []

        p_id = 0
        for obj in objects:
            points = self.clean_polygon(obj)
            shape = {
                'area': 0,
                'p_id': str(p_id),
                'points': [{'x': p[0], 'y': p[1]} for p in points]
            }

            area = nfp_utls.polygon_area(shape['points'])
            if area > 0:
                shape['points'].reverse()

            shape['area'] = abs(area)
            self.shapes.append(shape)

    def addContainer(self, container):

        """addContainer(object): adds a polygon objects as the container"""
        if not self.container:
            self.container = {}

        container = self.clean_polygon(container)

        self.container['points'] = [{'x': p[0], 'y':p[1]} for p in container]
        self.container['p_id'] = '-1'
        xbinmax = self.container['points'][0]['x']
        xbinmin = self.container['points'][0]['x']
        ybinmax = self.container['points'][0]['y']
        ybinmin = self.container['points'][0]['y']

        for point in self.container['points']:
            if point['x'] > xbinmax:
                xbinmax = point['x']
            elif point['x'] < xbinmin:
                xbinmin = point['x']
            if point['y'] > ybinmax:
                ybinmax = point['y']
            elif point['y'] < ybinmin:
                ybinmin = point['y']

        self.container['width'] = xbinmax - xbinmin
        self.container['height'] = ybinmax - ybinmin

    def clear(self):

        """clear(): Removes all objects and shape from the nester"""

        self.objects = None
        self.shapes = None

    def run(self):
        """run(): Runs a nesting operation. Returns a list of lists of
           shapes, each primary list being one filled container, or None
           if the operation failed."""
        # reset abort mechanism and variables

        # general conformity tests
        if not self.container:
            print("Empty container. Aborting")
            return
        if not self.shapes:
            print("Empty shapes. Aborting")
            return

        # and still identify the original face, so we can calculate a transform afterwards
        faces = list()
        for i in range(0, len(self.shapes)):
            shape = copy.deepcopy(self.shapes[i])
            shape['points'] = self.polygon_offset(shape['points'], self.config['spacing'])
            faces.append([
                str(i),
                shape
            ])
        # build a clean copy so we don't touch the original
        # order by area
        faces = sorted(faces, reverse=True, key=lambda face: face[1]['area'])
        self.bin_bounds = self.container['points']
        return self.launch_workers(faces)

    def launch_workers(self, adam):

        if self.GA is None:
            print('========== new ga')
            offset_bin = copy.deepcopy(self.container)
            offset_bin['points'] = self.polygon_offset(self.container['points'], self.config['spacing'])
            self.GA = geneticAlgorithm(adam, offset_bin, self.config)

        individual = None
        # 评估每个结果的适应值,多线程
        # for i in range(0, self.GA.config['populationSize']):
        #     if not self.GA.population[i].get('fitness'):
        #         individual = self.GA.population[i]
        #         print '===================GA population id %d' % i
        #     break
        #
        # if individual is None:
        #     # all individuals have been evaluated, start next generation
        #     self.GA.generation()
        #     individual = self.GA.population[1]
        res = list()
        for i in range(0, self.GA.config['populationSize']):
            res.append(self.find_fitness(self.GA.population[i]))
            self.GA.population[i]['fitness'] = res[-1]['fitness']
            # self.nfpCache.clear()

        self.show_result(res)

    def find_fitness(self, individual):
        place_list = individual['placement']
        rotations = individual['rotation']
        ids = [p[0] for p in place_list]
        for i in range(0, len(place_list)):
            place_list[i].append(rotations[i])

        nfp_pairs = list()
        newCache = dict()
        key = None
        for i in range(0, len(place_list)):
            part = place_list[i]
            key = {
                'A': '-1',
                'B': part[0],
                'inside': True,
                'A_rotation': 0,
                'B_rotation': rotations[i]
            }

            tmp_json_key = json.dumps(key)
            if not self.nfpCache.has_key(tmp_json_key):
                nfp_pairs.append({
                    'A': self.container,
                    'B': part[1],
                    'key': key
                })
            else:
                newCache[tmp_json_key] = self.nfpCache[tmp_json_key]

            # 没有算
            for j in range(0, i):
                placed = place_list[j]
                key = {
                    'A': placed[0],
                    'B': part[0],
                    'inside': False,
                    'A_rotation': rotations[j],
                    'B_rotation': rotations[i]
                }
                tmp_json_key = json.dumps(key)
                if not self.nfpCache.has_key(tmp_json_key):
                    nfp_pairs.append({
                        'A': placed[1],
                        'B': part[1],
                        'key': key
                    })
                else:
                    newCache[tmp_json_key] = self.nfpCache[tmp_json_key]

        # only keep cache for one cycle
        self.nfpCache = newCache
        # 另外的线程去排列图形
        if self.worker is None:
            self.worker = placement_worker.PlacementWorker(
                 self.container, place_list, ids, rotations, self.config, self.nfpCache)

        # TODO 多个基因的平行运算
        pair_list = list()
        for pair in nfp_pairs:
            pair_list.append(self.process_nfp(pair))

        # TODO:开子线程
        return self.generate_nfp(pair_list)

    def process_nfp(self, pair):
        if pair is None or len(pair) == 0:
            return None

        # 先考虑没有洞
        search_edges = self.config['exploreConcave']
        use_holes = self.config['useHoles']

        A = copy.deepcopy(pair['A'])
        A['points'] = nfp_utls.rotate_polygon(A['points'], pair['key']['A_rotation'])['points']
        B = copy.deepcopy(pair['B'])
        B['points'] = nfp_utls.rotate_polygon(B['points'], pair['key']['B_rotation'])['points']
        if pair['key']['inside']:
            if nfp_utls.is_rectangle(A['points'], 0.0001):
                nfp = nfp_utls.nfp_rectangle(A['points'], B['points'])
            else:
                nfp = nfp_utls.nfp_polygon(A, B, True, search_edges)

            # ensure all interior NFPs have the same winding direction
            if nfp and len(nfp) > 0:
                for i in range(0, len(nfp)):
                    if nfp_utls.polygon_area(nfp[i]) > 0:
                        nfp[i].reverse()
            else:
                print('NFP Warning:', pair['key'])

        else:
            if search_edges:
                nfp = nfp_utls.nfp_polygon(A, B, False, search_edges)
            else:
                nfp = minkowski_difference(A, B)

            # sanity check
            if nfp is None or len(nfp) == 0:
                print('error in NFP 260')
                # print('NFP Error:', pair['key'])
                # print('A;', A)
                # print('B:', B)
                return None

            for i in range(0, len(nfp)):
                # if searchedges is active, only the first NFP is guaranteed to pass sanity check
                if not search_edges or i == 0:
                    if abs(nfp_utls.polygon_area(nfp[i])) < abs(nfp_utls.polygon_area(A['points'])):
                        print('error in NFP area 269')
                        # print('NFP Area Error: ', abs(nfp_utls.polygon_area(nfp[i])), pair['key'])
                        # print('NFP:', json.dumps(nfp[i]))
                        # print('A: ', A)
                        # print('B: ', B)
                        nfp.pop(i)
                        return None

            if len(nfp) == 0:
                return None
            # for outer NFPs, the first is guaranteed to be the largest.
            # Any subsequent NFPs that lie inside the first are hole
            for i in range(0, len(nfp)):
                if nfp_utls.polygon_area(nfp[i]) > 0:
                    nfp[i].reverse()

                if i > 0:
                    if nfp_utls.point_in_polygon(nfp[i][0], nfp[0]):
                        if nfp_utls.polygon_area(nfp[i]) < 0:
                            nfp[i].reverse()

            # generate nfps for children (holes of parts) if any exist
            # TODO: 继续
            # 有洞的暂时不管
            if use_holes and len(A) > 0:
                pass
        return {'key': pair['key'], 'value': nfp}

    def generate_nfp(self, nfp):
        if nfp:
            for i in range(0, len(nfp)):
                if nfp[i]:
                    key = json.dumps(nfp[i]['key'])
                    self.nfpCache[key] = nfp[i]['value']

        self.worker.nfpCache.update(self.nfpCache)
        return self.worker.place_paths()

    def show_result(self, res):
        if len(res) > 0:
            best_result = res[0]

            for p in res:
                print p['fitness']
                if p['fitness'] < best_result['fitness']:
                    best_result = p

            if self.best is None or best_result['fitness'] > self.best['fitness']:
                self.best = best_result
                print self.best
                draw_svg(self.best['placements'], self.shapes, self.container)

    def polygon_offset(self, polygon, offset):
        is_list = True
        if isinstance(polygon[0], dict):
            polygon = [[p['x'], p['y']] for p in polygon]
            is_list = False

        miter_limit = 2
        co = pyclipper.PyclipperOffset(miter_limit, self.config['curveTolerance'])
        co.AddPath(polygon, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)
        result = co.Execute(1*offset)
        if not is_list:
            result = [{'x': p[0], 'y':p[1]} for p in result[0]]
        return result


    def clean_polygon(self, polygon):
        simple = pyclipper.SimplifyPolygon(polygon, pyclipper.PFT_NONZERO)

        if simple is None or len(simple) == 0:
            return None

        biggest = simple[0]
        biggest_area = pyclipper.Area(biggest)
        for i in range(1, len(simple)):
            area = abs(pyclipper.Area(simple[i]))
            if area > biggest_area:
                biggest = simple[i]
                biggest_area = area

        clean = pyclipper.CleanPolygon(biggest, self.config['curveTolerance'])
        if clean is None or len(clean) == 0:
            return None
        return clean


def draw_svg(shift_data, polygons, bin_polygon):
    if isinstance(polygons, Polygon):
        shapes = copy.deepcopy(polygons)
    else:
        shapes = list()
        for polygon in polygons:
            contour = [[p['x'], p['y']] for p in polygon['points']]
            shapes.append(Polygon(contour))

        shapes.append(Polygon([[p['x'], p['y']] for p in bin_polygon['points']]))

    # writeSVG('T3_2_0.svg', shapes, width=600)
    for s_data in shift_data:
        for move_step in s_data:
            if move_step['rotation'] <> 0:
                #center = shapes[int(move_step['p_id'])].center(0)
                shapes[int(move_step['p_id'])].rotate(math.pi / 180 * move_step['rotation'], 0, 0)
            shapes[int(move_step['p_id'])].shift(move_step['x'], move_step['y'])


    draw_polygon(shapes)
    #writeSVG('T3_6.svg', shapes, width=400)


class geneticAlgorithm():

    def __init__(self, adam, bin_polygon, config):
        self.bin_bounds = bin_polygon['points']
        self.bin_bounds = {
            'width': bin_polygon['width'],
            'height': bin_polygon['height'],
        }
        self.config = config
        self.bin_polygon = bin_polygon
        angles = list()
        for shape in adam:
            angles.append(self.random_angle(shape))

        self.population = [{'placement': adam, 'rotation': angles}]

        for i in range(1, self.config['populationSize']):
            mutant = self.mutate(self.population[0])
            self.population.append(mutant)

    def random_angle(self, shape):
        angle_list = list()
        for i in range(0, self.config['rotations']):
            angle_list.append(i * (360/self.config['rotations']))

        # 打乱顺序
        def shuffle_array(data):
            for i in range(len(data)-1, 0, -1):
                j = random.randint(0, i)
                data[i], data[j] = data[j], data[i]
            return data

        angle_list = shuffle_array(angle_list)

        # 查看选择后图形是否能放置在里面
        for angle in angle_list:
            rotate_part = nfp_utls.rotate_polygon(shape[1]['points'], angle)
            # 是否判断旋转出界,没有出界可以返回旋转角度,rotate 只是尝试去转，没有真正改变图形坐标
            if rotate_part['width'] < self.bin_bounds['width'] and rotate_part['height'] < self.bin_bounds['height']:
                return angle_list[i]

        return 0

    def mutate(self, individual):
        clone = {
            'placement': individual['placement'][:],
            'rotation': individual['rotation'][:]
        }
        for i in range(0, len(clone['placement'])):
            if random.random() < 0.01 * self.config['mutationRate']:
                if i+1 < len(clone['placement']):
                    clone['placement'][i],clone['placement'][i+1] = clone['placement'][i+1], clone['placement'][i]

        if random.random() < 0.01 * self.config['mutationRate']:
            clone['rotation'][i] = self.random_angle(clone['placement'][i])
        return clone

    def generation(self):
        # 适应度 从大到小排序
        self.population = sorted(self.population, reverse=True, key=lambda a: a['fitness'])
        new_population = [self.population[0]]
        while len(new_population) < self.config['populationSize']:
            male = self.random_weighted_individual()
            female = self.random_weighted_individual(male)
            # 交集下一代
            children = self.mate(male, female)

            # 轻微突变
            new_population.append(self.mutate(children[0]))

            if len(new_population) < self.config['populationSize']:
                new_population.append(self.mutate(children[1]))

        self.population = new_population

    def random_weighted_individual(self, exclude=None):
        pop = self.population
        if exclude and pop.index(exclude) >= 0:
            pop.remove(exclude)

        rand = random.random()
        lower = 0
        weight = 1 / len(pop)
        upper = weight
        pop_len = len(pop)
        for i in range(0, pop_len):
            if (rand > lower) and (rand < upper):
                return pop[i]
            lower = upper
            upper += 2 * weight * float(pop_len-i)/pop_len
        return pop[0]

    def mate(self, male, female):
        cutpoint = random.randint(0, len(male['placement'])-1)
        gene1 = male['placement'][:cutpoint]
        rot1 = male['rotation'][:cutpoint]

        gene2 = female['placement'][:cutpoint]
        rot2 = female['rotation'][:cutpoint]

        def contains(gene, shape_id):
            for i in range(0, len(gene)):
                if gene[i][0] == shape_id:
                    return True
            return False

        for i in range(len(female['placement']), 0, -1):
            if not contains(gene1, female['placement'][i][0]):
                gene1.append(female['placement'][i])
                rot1.append(female['rotation'][i])

        for i in range(len(male['placement']),0, -1):
            if not contains(gene2, male['placement'][i][0]):
                gene2.append(male['placement'][i])
                rot2.append(male['rotation'][i])

        return [{'placement': gene1, 'rotation': rot1}, {'placement': gene2, 'rotation': rot2}]


def minkowski_difference(A, B):
    Ac = [[p['x'], p['y']] for p in A['points']]
    Bc = [[p['x'] * -1, p['y'] * -1] for p in B['points']]
    solution = pyclipper.MinkowskiSum(Ac, Bc, True)
    largest_area = None
    clipper_nfp = None
    for p in solution:
        p = [{'x': i[0], 'y':i[1]} for i in p]
        sarea = nfp_utls.polygon_area(p)
        if largest_area is None or largest_area > sarea:
            clipper_nfp = p
            largest_area = sarea

    clipper_nfp = [{
                    'x': clipper_nfp[i]['x'] + Bc[0][0] * -1,
                    'y':clipper_nfp[i]['y'] + Bc[0][1] * -1
                   } for i in range(0, len(clipper_nfp))]
    return [clipper_nfp]


def draw_polygon(shapes):
    ax = plt.axes()
    ax.set_xlim(-10, 1500)
    ax.set_ylim(-10, 1500)
    output_obj = list()
    output_obj.append(patches.Polygon(shapes[-1].contour(0), fc='green'))
    print shapes[-1].contour(0)
    shapes.pop(-1)

    for s in shapes:
        print s.contour(0)
        # x, y ,w, h =  s.boundingBox()
        # output_obj.append(patches.Rectangle((x, y), w, h, fc='green'))
        output_obj.append(patches.Polygon(s.contour(0), fc='yellow', lw=1, edgecolor='m'))
    # corner_points_2, hull_points 是最终坐标，然后通过新的坐标
    for p in output_obj:
        ax.add_patch(p)
    plt.gca().set_aspect('equal')
    plt.title('Rate')
    plt.show()

if __name__ == '__main__':
    n = Nester()
    s = input_utls.input_polygon('dxf_file/T4.dxf', is_class=False)
    # writeSVG('T2_O.svg', s, width=600)

    n.addObjects(s[1:])
    n.addContainer(s[0])
    n.run()
    # n.show()
