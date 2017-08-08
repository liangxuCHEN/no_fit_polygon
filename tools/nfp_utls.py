# -*- coding: utf-8 -*-
import copy
import math
TOL = 0.0000001   # 计算过程中误差忽略值


def almost_equal(a, b, tolerance=None):
    """
    A，B 两点是否约为相同
    :param a: 坐标
    :param b: 坐标
    :param tolerance: 误差值
    :return:
    """
    if tolerance is None:
        tolerance = TOL
    return abs(a - b) < tolerance


def is_rectangle(poly, tolerance=None):
    bb = get_polygon_bounds(poly)
    tolerance = tolerance or TOL
    for point in poly:
        if not almost_equal(point['x'], bb['x'], tolerance) and not almost_equal(
                point['x'], bb['x'] + bb['width'], tolerance):
            return False
        if not almost_equal(point['y'], bb['y'], tolerance) and not almost_equal(
                point['y'], bb['y'] + bb['height'], tolerance):
            return False

    return True


def normalize_vector(v):
    """
    normalize vector into a unit vector
    :return:
    """
    if almost_equal(v['x'] * v['x'] + v['y'] * v['y'], 1):
        # given vector was already a unit vector
        return v
    inverse = 1.0 / math.sqrt(v['x']**2 + v['y']**2)

    return {'x': v['x']*inverse, 'y': v['y']*inverse}


def on_segment(A, B, p):
    """
    returns true if p lies on the line segment defined by AB, but not at any endpoints
    :param A:
    :param B:
    :param p:
    :return:
    """
    # vertical line
    if almost_equal(A['x'], B['x']) and almost_equal(p['x'], A['x']):
        if not almost_equal(p['y'], B['y']) and not almost_equal(p['y'], A['y']) and \
                        max(B['y'], A['y']) > p['y'] and p['y'] > min(B['y'], A['y']):
            return True
        else:
            return False
    # vertical line
    if almost_equal(A['y'], B['y']) and almost_equal(p['y'], A['y']):
        if not almost_equal(p['x'], B['x']) and not almost_equal(p['x'], A['x']) and \
                        max(B['x'], A['x']) > p['x'] and p['x'] > min(B['x'], A['x']):
            return True
        else:
            return False
    # range check
    if (p['x'] < A['x'] and p['x'] < B['x']) or (p['x'] > A['x'] and p['x'] > B['x']) or (
                    p['y'] < A['y'] and p['y'] < B['y']) or (p['y'] > A['y'] and p['y'] > B['y']):
        return False

    # exclude end points
    if (almost_equal(p['x'], A['x']) and almost_equal(p['y'], A['y'])) or (
                almost_equal(p['x'], B['x']) and almost_equal(p['y'], B['y'])):
        return False

    cross = (p['y'] - A['y']) * (B['x'] - A['x']) - (p['x'] - A['x']) * (B['y'] - A['y'])
    if abs(cross) > TOL:
        return False
    dot = (p['x'] - A['x']) * (B['x'] - A['x']) + (p['y'] - A['y']) * (B['y'] - A['y'])
    if dot < 0 or almost_equal(dot, 0):
        return False

    len2 = (B['x'] - A['x']) * (B['x'] - A['x']) + (B['y'] - A['y']) * (B['y'] - A['y'])
    if dot > len2 or almost_equal(dot, len2):
        return False
    return True


def nfp_rectangle(A, B):
    """
    :param A: {x:12, y:10}
    :param B: {x:12, y:10}
    :return:
    """
    min_ax = A[0]['x']
    min_ay = A[0]['y']
    max_ax = A[0]['x']
    max_ay = A[0]['y']

    for point in A[1:]:
        if point['x'] < min_ax:
            min_ax = point['x']
        if point['x'] > max_ax:
            max_ax = point['x']
        if point['y'] < min_ay:
            min_ay = point['y']
        if point['y'] > max_ay:
            max_ay = point['y']

    min_bx = B[0]['x']
    min_by = B[0]['y']
    max_bx = B[0]['x']
    max_by = B[0]['y']

    for point in B[1:]:
        if point['x'] < min_bx:
            min_bx = point['x']
        if point['x'] > max_bx:
            max_bx = point['x']
        if point['y'] < min_by:
            min_by = point['y']
        if point['y'] > max_by:
            max_by = point['y']

    if max_bx - min_bx > max_ax - min_ax:
        return None
    if max_by - min_by > max_ay - min_ay:
        return None

    return [[
        {'x': min_ax-min_bx+B[0]['x'], 'y': min_ay-min_by+B[0]['y']},
        {'x': max_ax-max_bx+B[0]['x'], 'y': min_ay-min_by+B[0]['y']},
        {'x': max_ax-max_bx+B[0]['x'], 'y': max_ay-max_by+B[0]['y']},
        {'x': min_ax-min_bx+B[0]['x'], 'y': max_ay-max_by+B[0]['y']}
    ]]


def nfp_polygon(A, B, inside=True, search_edges=False):
    """
    given a static polygon A and a movable polygon B, compute a no fit polygon by orbiting B about A
    if the inside flag is set, B is orbited inside of A rather than outside
    if the searchEdges flag is set, all edges of A are explored for NFPs - multiple
    """
    if A is None or len(A['points']) < 3 or B is None or len(B['points']) < 3:
        return None

    # A last point = offsetx, offsety
    A['offsetx'] = 0
    A['offsety'] = 0

    min_a = A['points'][0]['y']
    min_a_index = 0
    max_b = B['points'][0]['y']
    max_b_index = 0

    for i in range(1, len(A['points'])):
        A['points'][i]['marked'] = False
        if A['points'][i]['y'] < min_a:
            min_a = A['points'][i]['y']
            min_a_index = i

    for i in range(1, len(B['points'])):
        B['points'][i]['marked'] = False
        if B['points'][i]['y'] > max_b:
            max_b = B['points'][i]['y']
            max_b_index = i

    if not inside:
        # shift B such that the bottom-most point of B is at the top-most point of A.
        # This guarantees an initial placement with no intersections
        start_point = {
            'x': A['points'][min_a_index]['x'] - B['points'][max_b_index]['x'],
            'y': A['points'][min_a_index]['y'] - B['points'][max_b_index]['y']
        }
    else:
        #  no reliable heuristic for inside
        start_point = search_start_point(A, B, inside)

    NFP_list = list()

    while start_point:
        B['offsetx'] = start_point['x']
        B['offsety'] = start_point['y']

        # maintain a list of touching points/edges
        prevvector = None
        NFP = [{
            'x': B['points'][0]['x'] + B['offsetx'],
            'y': B['points'][0]['y'] + B['offsety'],
        }]

        referencex = B['points'][0]['x'] + B['offsetx']
        referencey = B['points'][0]['y'] + B['offsety']
        startx = referencex
        starty = referencey
        counter = 0
        len_a = len(A['points'])
        len_b = len(B['points'])
        while counter < 10 * (len_a + len_b):
            touching = list()
            for i in range(0, len_a):
                nexti = 0 if i == len_a - 1 else i + 1
                for j in range(len_b):
                    nextj = 0 if j == len_b - 1 else j + 1
                    if almost_equal(A['points'][i]['x'], B['points'][j]['x'] + B['offsetx']) and almost_equal(
                            A['points'][i]['y'], B['points'][j]['y'] + B['offsety']):
                        touching.append({'type': 0, 'A': i, 'B': j})
                    elif on_segment(A['points'][i], A['points'][nexti],
                                    {'x': B['points'][j]['x']+B['offsetx'], 'y': B['points'][j]['y'] + B['offsety']}):
                        touching.append({'type': 1, 'A': nexti, 'B': j})
                    elif on_segment(
                            {'x': B['points'][j]['x']+B['offsetx'], 'y': B['points'][j]['y'] + B['offsety']},
                            {'x': B['points'][nextj]['x'] + B['offsetx'], 'y': B['points'][nextj]['y'] + B['offsety']},
                            A['points'][i]):
                        touching.append({'type': 2, 'A': i, 'B': nextj})

            # generate translation vectors from touching vertices/edges
            vectors = list()
            for i in range(0, len(touching)):
                vertex_a = {'A': A['points'][touching[i]['A']], 'marked': True}
                prev_a_index = touching[i]['A'] - 1
                next_a_index = touching[i]['A'] + 1
                prev_a_index = len_a - 1 if prev_a_index < 0 else prev_a_index  # loop
                next_a_index = 0 if next_a_index >= len_a else next_a_index  # loop

                prev_a = A['points'][prev_a_index]
                next_a = A['points'][next_a_index]

                # adjacent B vertices
                vertex_b = {'A': B['points'][touching[i]['B']]}
                prev_b_index = touching[i]['B'] - 1
                next_b_index = touching[i]['B'] + 1
                prev_b_index = len_b - 1 if prev_b_index < 0 else prev_b_index  # loop
                next_b_index = 0 if next_b_index >= len_b else next_b_index  # loop

                prev_b = B['points'][prev_b_index]
                next_b = B['points'][next_b_index]

                if touching[i]['type'] == 0:
                    v_a1 = {
                        'x': prev_a['x'] - vertex_a['A']['x'],
                        'y': prev_a['y'] - vertex_a['A']['y'],
                        'start': vertex_a['A'],
                        'end': prev_a
                    }

                    v_a2 = {
                        'x': next_a['x'] - vertex_a['A']['x'],
                        'y': next_a['y'] - vertex_a['A']['y'],
                        'start': vertex_a['A'],
                        'end': next_a
                    }

                    v_b1 = {
                        'x': vertex_b['A']['x'] - prev_b['x'],
                        'y': vertex_b['A']['y'] - prev_b['y'],
                        'start': prev_b,
                        'end': vertex_b['A']
                    }

                    v_b2 = {
                        'x': vertex_b['A']['x'] - next_b['x'],
                        'y': vertex_b['A']['y'] - next_b['y'],
                        'start': next_b,
                        'end': vertex_b['A']
                    }

                    vectors.append(v_a1)
                    vectors.append(v_a2)
                    vectors.append(v_b1)
                    vectors.append(v_b2)
                elif touching[i]['type'] == 1:
                    vectors.append({
                        'x': vertex_a['A']['x'] - (vertex_b['A']['x'] + B['offsetx']),
                        'y': vertex_a['A']['y'] - (vertex_b['A']['y'] + B['offsety']),
                        'start': prev_a,
                        'end': vertex_a['A']
                    })
                    vectors.append({
                        'x': prev_a['x'] - (vertex_b['A']['x'] + B['offsetx']),
                        'y': prev_a['y'] - (vertex_b['A']['y'] + B['offsety']),
                        'start': vertex_a['A'],
                        'end': prev_a
                    })
                elif touching[i]['type'] == 2:
                    vectors.append({
                        'x': vertex_a['A']['x'] - (vertex_b['A']['x'] + B['offsetx']),
                        'y': vertex_a['A']['y'] - (vertex_b['A']['y'] + B['offsety']),
                        'start': prev_b,
                        'end': vertex_b['A']
                    })
                    vectors.append({
                        'x': vertex_a['A']['x'] - (prev_b['x'] + B['offsetx']),
                        'y': vertex_a['A']['y'] - (prev_b['y'] + B['offsety']),
                        'start': vertex_b['A'],
                        'end': prev_b
                    })

            translate = None
            max_d = 0
            for i in range(0, len(vectors)):
                if vectors[i]['x'] == 0 and vectors[i]['y'] == 0:
                    continue
                # if this vector points us back to where we came from, ignore it.
                # ie cross product = 0, dot product < 0
                if prevvector and (vectors[i]['y'] * prevvector['y'] + vectors[i]['x'] * prevvector['x']) < 0:
                    # compare magnitude with unit vectors
                    vectorlength = math.sqrt(vectors[i]['x']**2 + vectors[i]['y']**2)
                    unitv = {'x': vectors[i]['x'] / vectorlength, 'y': vectors[i]['y'] / vectorlength}
                    prevlength = math.sqrt(prevvector['x']**2+prevvector['y']**2)
                    prevunit = {'x': prevvector['x'] / prevlength, 'y': prevvector['y'] / prevlength}

                    # we need to scale down to unit vectors to normalize vector length. Could also just do a tan here
                    if abs(unitv['y'] * prevunit['x'] - unitv['x'] * prevunit['y']) < 0.0001:
                        continue

                d = polygon_slide_distance(A, B, vectors[i], True)
                vecd2 = vectors[i]['x']**2 + vectors[i]['y']**2

                if d is None or d**2 > vecd2:
                    vecd = math.sqrt(vectors[i]['x']**2 + vectors[i]['y']**2)
                    d = vecd
                if d and d > max_d:
                    max_d = d
                    translate = vectors[i]

            if translate is None or almost_equal(max_d, 0):
                NFP = None
                break

            translate['start']['marked'] = True
            translate['end']['marked'] = True

            prevvector = translate

            # trim
            vlength2 = translate['x']**2 + translate['y']**2
            if max_d**2 < vlength2 and not almost_equal(max_d**2, vlength2):
                scale = math.sqrt((max_d**2)/vlength2)
                translate['x'] *= scale
                translate['y'] *= scale

            referencex += translate['x']
            referencey += translate['y']

            if almost_equal(referencex, startx) and almost_equal(referencey, starty):
                # we have made a full loop
                break

            # if A and B start on a touching horizontal line, the end point may not be the start point
            looped = False
            if len(NFP) > 0:
                for i in range(0, len(NFP)-1):
                    if almost_equal(referencex, NFP[i]['x'] and almost_equal(referencey, NFP[i]['y'])):
                        looped = True

            if looped:
                # we've made a full loop
                break

            NFP.append({
                'x': referencex,
                'y': referencey
            })
            B['offsetx'] += translate['x']
            B['offsety'] += translate['y']

            counter += 1

            if NFP and len(NFP) > 0:
                NFP_list.append(NFP)

        if not search_edges:
            # only get outer NFP or first inner NFP
            break

        start_point = search_start_point(A, B, inside, NFP_list)

    return NFP_list


def search_start_point(A, B, inside=True, NFP=None):
    """
    searches for an arrangement of A and B such that they do not overlap if an NFP is given,
    only search for startpoints that have not already been traversed in the given NFP
    :param A:
    :param B:
    :param inside:
    :param NFP:
    :return:
    """
    # clone arrays
    A = copy.deepcopy(A)
    B = copy.deepcopy(B)

    for i in range(0, len(A['points'])-1):
        if not A['points'][i].get('marked'):
            A['points'][i]['marked'] = True
            for j in range(0, len(B['points'])):
                B['offsetx'] = A['points'][i]['x'] - B['points'][j]['x']
                B['offsety'] = A['points'][i]['y'] - B['points'][j]['y']

                # 判断 A，B 是否一样
                # 点是否在多边形边上
                bin_side = None
                for k in range(0, len(B['points'])):
                    inpoly = point_in_polygon(
                        {'x': B['points'][k]['x']+B['offsetx'],
                         'y': B['points'][k]['y']+B['offsety']}, A)
                    if inpoly is not None:
                        bin_side = inpoly
                        break

                if bin_side is None:
                    return None

                start_point = {
                    'x': B['offsetx'],
                    'y': B['offsety']
                }
                if ((bin_side and inside) or (not bin_side and not inside)) and (
                            not intersect(A, B) and not inNfp(start_point, NFP)):
                    return start_point

                # slide B along vector
                vx = A['points'][i+1]['x'] - A['points'][i]['x']
                vy = A['points'][i+1]['y'] - A['points'][i]['y']

                d1 = polygon_projection_distance(A, B, {'x': vx, 'y': vy})
                d2 = polygon_projection_distance(B, A, {'x': -vx, 'y': -vy})

                d = None

                if d1 is not None and d2 is not None:
                    d = min(d1, d2)
                elif d1 is None and d2 is not None:
                    d = d2
                elif d1 is not None and d2 is None:
                    d = d1

                # only slide until no longer negative
                if not (d is not None and not almost_equal(d, 0) and d > 0):
                    continue

                vd2 = vx * vx + vy * vy
                if d * d < vd2 and not almost_equal(d*d, vd2):
                    vd = math.sqrt(vx * vx + vy * vy)
                    vx *= d /vd
                    vy *= d /vd

                B['offsetx'] += vx
                B['offsety'] += vy

                for k in range(0, len(B['points'])):
                    inpoly = point_in_polygon(
                        {'x': B['points'][k]['x']+B['offsetx'],
                         'y': B['points'][k]['y']+B['offsety']}, A)
                    if inpoly is not None:
                        bin_side = inpoly
                        break

                start_point = {'x': B['offsetx'], 'y': B['offsety']}
                if ((bin_side and inside) or (not bin_side and not inside)) and (
                            not intersect(A, B) and not inNfp(start_point, NFP)):
                    return start_point

    return None


def inNfp(p, nfp):
    """
    returns true if point already exists in the given nfp
    :param p:
    :param nfp:
    :return:
    """
    if not nfp or len(nfp) == 0:
        return False

    for i in range(0, len(nfp)):
        for j in range(0, len(nfp[i])):
            if almost_equal(p['x'], nfp[i][j]['x']) and almost_equal(p['y'], nfp[i][j]['y']):
                return True

    return False


def point_in_polygon(point, polygon):
    if isinstance(polygon, list):
        polygon = {'points': polygon}
    if len(polygon.get('points')) < 3:
        return None

    inside = False
    offsetx = polygon.get('offsetx') or 0
    offsety = polygon.get('offsety') or 0

    j = len(polygon['points']) - 1
    for i in range(0, len(polygon['points'])):
        xi = polygon['points'][i]['x'] + offsetx
        yi = polygon['points'][i]['y'] + offsety
        xj = polygon['points'][j]['x'] + offsetx
        yj = polygon['points'][j]['y'] + offsety

        if almost_equal(xi, point['x']) and almost_equal(yi, point['y']):
            return None

        if on_segment({'x': xi, 'y': yi}, {'x':xj, 'y':yj}, point):
            return None  # exactly on the segment

        if almost_equal(xi, xj) and almost_equal(yi, yj):
            # ignore very small lines
            continue

        intersect = ((yi > point['y']) != (yj > point['y'])) and (point['x'] < (xj - xi) * (point['y'] - yi) / (yj - yi) + xi)
        if intersect:
            inside = not inside

    return inside


def intersect(A, B):

    a_offsetx = A['offsetx'] or 0
    a_offsety = A['offsety'] or 0

    b_offsetx = B['offsetx'] or 0
    b_offsety = B['offsety'] or 0

    A = copy.deepcopy(A)
    B = copy.deepcopy(B)
    len_a = len(A['points'])
    len_b = len(B['points'])
    for i in range(0, len_a - 1):
        for j in range(0, len_b - 1):
            a1 = {'x': A['points'][i]['x']+a_offsetx, 'y': A['points'][i]['y']+a_offsety}
            a2 = {'x': A['points'][i+1]['x']+a_offsetx, 'y': A['points'][i+1]['y']+a_offsety}
            b1 = {'x': B['points'][j]['x']+b_offsetx, 'y': B['points'][j]['y']+b_offsety}
            b2 = {'x': B['points'][j+1]['x']+b_offsetx, 'y': B['points'][j+1]['y']+b_offsety}

            pre_vb_index = len_b - 1 if j == 0 else j - 1
            pre_va_index = len_a - 1 if i == 0 else i - 1
            next_b_index = 0 if j + 1 == len_b - 1 else j + 2
            next_a_index = 0 if i + 1 == len_a - 1 else i + 2

            # go even further back if we happen to hit on a loop end point
            if B['points'][pre_vb_index] == B['points'][j] or almost_equal(
                    B['points'][pre_vb_index]['x'], B['points'][j]['x']) and almost_equal(
                B['points'][pre_vb_index]['y'], B['points'][j]['y']):
                pre_vb_index = len_b - 1 if pre_vb_index == 0 else pre_vb_index - 1

            if A['points'][pre_va_index] == A['points'][i] or almost_equal(
                    A['points'][pre_va_index]['x'], A['points'][i]['x']) and almost_equal(
                A['points'][pre_va_index]['y'], A['points'][i]['y']):
                pre_va_index = len_a - 1 if pre_va_index == 0 else pre_va_index - 1

            # go even further forward if we happen to hit on a loop end point
            if B['points'][next_b_index] == B['points'][j+1] or almost_equal(
                    B['points'][next_b_index]['x'], B['points'][j+1]['x']) and almost_equal(
                B['points'][next_b_index]['y'], B['points'][j+1]['y']):
                next_b_index = 0 if next_b_index == len_b - 1 else next_b_index + 1

            if A['points'][next_a_index] == A['points'][i+1] or almost_equal(
                    A['points'][next_a_index]['x'], A['points'][i+1]['x']) and almost_equal(
                A['points'][next_a_index]['y'], A['points'][i+1]['y']):
                next_a_index = 0 if next_a_index == len_a - 1 else next_a_index + 1

            a0 = {'x': A['points'][pre_va_index]['x']+a_offsetx, 'y': A['points'][pre_va_index]['y']+a_offsety}
            b0 = {'x': B['points'][pre_vb_index]['x']+b_offsetx, 'y': B['points'][pre_vb_index]['y']+b_offsety}
            a3 = {'x': A['points'][next_a_index]['x']+a_offsetx, 'y': A['points'][next_a_index]['y']+a_offsety}
            b3 = {'x': B['points'][next_b_index]['x']+b_offsetx, 'y': B['points'][next_b_index]['y']+b_offsety}

            if on_segment(a1, a2, b1) or almost_equal(a1['x'], b1['x']) and almost_equal(a1['y'], b1['y']):
                # if a point is on a segment, it could intersect or it could not. Check via the neighboring points
                b0in = point_in_polygon(b0, A)
                b2in = point_in_polygon(b2, A)
                if (b0in and not b2in) or (not b0in and b2in):
                    return True
                else:
                    continue

            if on_segment(a1, a2, b2) or almost_equal(a2['x'], b2['x']) and almost_equal(a2['y'], b2['y']):
                # if a point is on a segment, it could intersect or it could not.Check via the neighboring points
                b1in = point_in_polygon(b1, A)
                b3in = point_in_polygon(b3, A)
                if (b1in and not b3in) or (not b1in and b3in):
                    return True
                else:
                    continue

            if on_segment(b1, b2, a1) or almost_equal(a1['x'], b2['x']) and almost_equal(a1['y'], b2['y']):
                # if a point is on a segment, it could intersect or it could not.Check via the neighboring points
                a0in = point_in_polygon(a0, B)
                a2in = point_in_polygon(a2, B)
                if (a0in and not a2in) or (not a0in and a2in):
                    return True
                else:
                    continue

            if on_segment(b1, b2, a2) or almost_equal(a2['x'], b1['x']) and almost_equal(a2['y'], b1['y']):
                # if a point is on a segment, it could intersect or it could not.Check via the neighboring points
                a1in = point_in_polygon(a1, B)
                a3in = point_in_polygon(a3, B)
                if (a1in and not a3in) or (not a1in and a3in):
                    return True
                else:
                    continue

            if line_intersect(b1, b2, a1, a2):
                return True

    return False


def line_intersect(A, B, E, F, infinite=None):
    """
    returns the intersection of AB and EF, or null if there are no intersections or other numerical error
    if the infinite flag is set, AE and EF describe infinite lines without endpoints,
    they are finite line segments otherwise
    :param A:
    :param B:
    :param E:
    :param F:
    :param infinite:
    :return:
    """
    a1 = B['y'] - A['y']
    b1 = A['x'] - B['x']
    c1 = B['x'] * A['y'] - A['x'] * B['y']
    a2 = F['y'] - E['y']
    b2 = E['x'] - F['y']
    c2 = F['x'] * E['y'] - E['x'] * F['y']
    denom = a1 * b2 - a2 * b1
    if denom == 0:
        return None
    x = (b1 * c2 - b2 * c1) / denom
    y = (a2 * c1 - a1 * c2) / denom

    if infinite is None:
        if abs(A['x'] - B['x']) > TOL:
            tmp = x < A['x'] or x > B['x'] if A['x'] < B['x'] else x > A['x'] or x < B['x']
            if tmp:
                return None
            tmp = y < A['y'] or y > B['y'] if A['y'] < B['y'] else y > A['y'] or y < B['y']
            if tmp:
                return None
        if abs(E['x'] - F['x']) > TOL:
            tmp = x < E['x'] or x > F['x'] if E['x'] < F['x'] else x > E['x'] or x < F['x']
            if tmp:
                return None
            tmp = y < E['y'] or y > F['y'] if E['y'] < F['y'] else y > E['y'] or y < F['y']
            if tmp:
                return None

    return {'x': x, 'y': y}


def polygon_projection_distance(A, B, direction):
    """
    project each point of B onto A in the given direction, and return the distance
    :param A:
    :param B:
    :param direction:
    :return:
    """
    b_offsetx = B.get('offsetx') or 0
    b_offsety = B.get('offsety') or 0
    a_offsetx = A.get('offsetx') or 0
    a_offsety = A.get('offsety') or 0

    A = copy.deepcopy(A)
    B = copy.deepcopy(B)

    edge_a = A['points']
    edge_b = B['points']
    distance = None
    p = dict()
    s1 = dict()
    s2 = dict()
    for i in range(0, len(edge_b)):
        # the shortest/most negative projection of B onto A
        min_projection = minp = None
        for j in range(0, len(edge_a) - 1):
            p['x'] = edge_b[i]['x'] + b_offsetx
            p['y'] = edge_b[i]['y'] + b_offsety
            s1['x'] = edge_a[j]['x'] + a_offsetx
            s1['y'] = edge_a[j]['y'] + a_offsety
            s2['x'] = edge_a[j+1]['x'] + a_offsetx
            s2['y'] = edge_a[j+1]['y'] + a_offsety

            if abs((s2['y'] - s1['y']) * direction['x'] - (s2['x'] - s2['x']) * direction['y']) < TOL:
                continue

            # project point, ignore edge boundaries
            d = point_distance(p, s1, s2, direction)
            if d and (min_projection is None or d < min_projection):
                min_projection = d

        if min_projection and (distance is None or min_projection > distance):
            distance = min_projection

    return distance


def point_distance(p, s1, s2, normal, infinite=None):
    normal = normalize_vector(normal)

    dir_point = {
        'x': normal['y'],
        'y': -normal['x'],
    }

    pdot = p['x'] * dir_point['x'] + p['y'] * dir_point['y']
    s1dot = s1['x'] * dir_point['x'] + s1['y'] * dir_point['y']
    s2dot = s2['x'] * dir_point['x'] + s2['y'] * dir_point['y']

    pdotnorm = p['x']*normal['x'] + p['y'] * normal['y']
    s1dotnorm = s1['x']*normal['x'] + s1['y'] * normal['y']
    s2dotnorm = s2['x'] * normal['x'] + s2['y'] * normal['y']

    if infinite is None:
        # dot doesn't collide with segment, or lies directly on the vertex
        if ((pdot < s1dot or almost_equal(pdot, s1dot)) and (pdot < s2dot or almost_equal(pdot, s2dot))) or (
                    (pdot > s1dot or almost_equal(pdot, s1dot)) and ((pdot > s2dot) or almost_equal(pdot, s2dot))):
            return None
        if (almost_equal(pdot, s1dot) and almost_equal(pdot, s2dot)) and (
                        pdotnorm > s1dotnorm and pdotnorm > s2dotnorm):
            return min(pdotnorm - s1dotnorm, pdotnorm - s2dotnorm)
        if almost_equal(pdot, s1dot) and almost_equal(pdot, s2dot) and pdotnorm < s1dotnorm and pdotnorm < s2dotnorm:
            return -min(s1dotnorm-pdotnorm, s2dotnorm-pdotnorm)

    return -(pdotnorm - s1dotnorm + (s1dotnorm - s2dotnorm) * (s1dot - pdot)/(s1dot - s2dot))


def polygon_slide_distance(A, B, direction, ignorenegative):

    b_offsetx = B.get('offsetx') or 0
    b_offsety = B.get('offsety') or 0
    a_offsetx = A.get('offsetx') or 0
    a_offsety = A.get('offsety') or 0

    A = copy.deepcopy(A)
    B = copy.deepcopy(B)

    if not A['points'][-1] == A['points'][0]:
        A['points'].append(A['points'][0])
    if not B['points'][0] == B['points'][-1]:
        B['points'].append(B['points'][0])

    edge_a = A['points']
    edge_b = B['points']
    distance = None

    dir_point = normalize_vector(direction)

    for i in range(0, len(edge_b) - 1):

        for j in range(0, len(edge_a) - 1):
            A1 = {'x': edge_a[j]['x'] + a_offsetx, 'y': edge_a[j]['y'] + a_offsety}
            A2 = {'x': edge_a[j+1]['x'] + a_offsetx, 'y': edge_a[j+1]['y'] + a_offsety}
            B1 = {'x': edge_b[i]['x'] + b_offsetx, 'y': edge_b[i]['y'] + b_offsety}
            B2 = {'x': edge_b[i + 1]['x'] + b_offsetx, 'y': edge_b[i + 1]['y'] + b_offsety}

            if (almost_equal(A1['x'], A2['x']) and almost_equal(A1['y'], A2['y'])) or almost_equal(
                    B1['x'], B2['x']) and almost_equal(B1['y'], B2['y']):
                continue

            d = segment_distance(A1, A2, B1, B2, dir_point)
            if d and (distance is None or d < distance):
                if not ignorenegative or d > 0 or almost_equal(d, 0):
                    distance = d
    return distance


def segment_distance(A, B, E, F, direction):
    normal = {
        'x': direction['y'],
        'y': -direction['x']
    }
    reverse = {
        'x': -direction['x'],
        'y': -direction['y']
    }

    dot_a = A['x'] * normal['x'] + A['y'] * normal['y']
    dot_b = B['x'] * normal['x'] + B['y'] * normal['y']
    dot_e = E['x'] * normal['x'] + E['y'] * normal['y']
    dot_f = F['x'] * normal['x'] + F['y'] * normal['y']

    cross_a = A['x'] * direction['x'] + A['y'] * direction['y']
    cross_b = B['x'] * direction['x'] + B['y'] * direction['y']
    cross_e = E['x'] * direction['x'] + E['y'] * direction['y']
    cross_f = F['x'] * direction['x'] + F['y'] * direction['y']

    ab_min = min(dot_a, dot_b)
    ab_max = max(dot_a, dot_b)

    ef_min = min(dot_e, dot_f)
    ef_max = max(dot_e, dot_f)

    # segments that will touch at one point
    if almost_equal(ab_max, ef_min, TOL) or almost_equal(ab_min, ef_max, TOL):
        return None

    # segments miss each other completely
    if ab_max < ef_min or ab_min > ef_max:
        return None

    if (ab_max > ef_max and ab_min < ef_min) or (ef_max > ab_max and ef_min < ab_min):
        overlap = 1
    else:
        min_max = min(ab_max, ef_max)
        max_min = max(ab_min, ef_min)
        max_max = max(ab_max, ef_max)
        min_min = min(ab_min, ef_min)

        overlap = (min_max - max_min) / (max_max - min_min)

    cross_abe = (E['y'] - A['y']) * (B['x'] - A['x']) - (E['x'] - A['x']) * (B['y'] - A['y'])
    cross_abf = (F['y'] - A['y']) * (B['x'] - A['x']) - (F['x'] - A['x']) * (B['y'] - A['y'])

    # lines are colinear
    if almost_equal(cross_abe, 0) and almost_equal(cross_abf, 0):
        ab_norm = {'x': B['y'] - A['y'], 'y': A['x'] - B['x']}
        ef_norm = {'x': F['y'] - E['y'], 'y': E['x'] - F['x']}

        ab_norm_length = math.sqrt(ab_norm['x']**2 + ab_norm['y']**2)
        ab_norm['x'] /= ab_norm_length
        ab_norm['y'] /= ab_norm_length

        ef_norm_length = math.sqrt(ef_norm['x']**2 + ef_norm['y']**2)
        ef_norm['x'] /= ef_norm_length
        ef_norm['y'] /= ef_norm_length

        # segment normals must point in opposite directions
        if abs(ab_norm['y'] * ef_norm['x'] - ab_norm['x'] * ef_norm['y']) < TOL and (
                            ab_norm['y'] * ef_norm['y'] + ab_norm['x'] * ef_norm['x'] < 0):
            # normal of AB segment must point in same direction as given direction vector
            norm_dot = ab_norm['y'] * direction['y'] + ab_norm['x'] * direction['x']
            # the segments merely slide along eachother
            if almost_equal(norm_dot, 0, TOL):
                return None

            if norm_dot < 0:
                return 0

        return None

    distances = list()

    # coincident points
    if almost_equal(dot_a, dot_e):
        distances.append(cross_a - cross_e)
    elif almost_equal(dot_a, dot_f):
        distances.append(cross_a - cross_f)
    elif ef_min < dot_a and dot_a < ef_max:
        d = point_distance(A, E, F, reverse)
        # A currently touches EF, but AB is moving away from EF
        if d and almost_equal(d, 0):
            db = point_distance(B, E, F, reverse, True)
            if db < 0 or almost_equal(db * overlap, 0):
                d = None
        if d:
            distances.append(d)

    if almost_equal(dot_b, dot_e):
        distances.append(cross_b - cross_e)
    elif almost_equal(dot_b, dot_f):
        distances.append(cross_b - cross_f)
    elif dot_b > ef_min and dot_b < ef_max:
        d = point_distance(B, E, F, reverse)

        if d and almost_equal(d, 0):
            da = point_distance(A, E, F, reverse, True)
            if da < 0 or almost_equal(da * overlap, 0):
                d = None

        if d:
            distances.append(d)

    if dot_e > ab_min and dot_e < ab_max:
        d = point_distance(E, A, B, direction)
        if d and almost_equal(d, 0):
            df = point_distance(F, A, B, direction, True)
            if df < 0 or almost_equal(df * overlap, 0):
                d = None
        if d:
            distances.append(d)

    if len(distances) == 0:
        return None

    if dot_f > ab_min and dot_f < ab_max:
        d = point_distance(F, A, B, direction)
        if d and almost_equal(d, 0):
            de = point_distance(E, A, B, direction, True)
            if de < 0 or almost_equal(de * overlap, 0):
                d = None
        if d:
            distances.append(d)

    if len(distances) == 0:
        return None

    return min(distances)


def polygon_area(polygon):
    area = 0
    j = len(polygon) - 1
    for i in range(0, len(polygon)):
        area += (polygon[j]['x'] + polygon[i]['x']) * (polygon[j]['y'] - polygon[i]['y'])
        j = i

    return 0.5 * area


def rotate_polygon(polygon, angle):
    rotated = {'points': list()}
    angle = angle * math.pi / 180
    for p in polygon:
        x = p['x']
        y = p['y']
        rotated['points'].append({
            'x': x * math.cos(angle) - y * math.sin(angle),
            'y': x * math.sin(angle) + y * math.cos(angle)
        })

    bounds = get_polygon_bounds(rotated['points'])
    rotated['x'] = bounds['x']
    rotated['y'] = bounds['y']
    rotated['width'] = bounds['width']
    rotated['height'] = bounds['height']
    return rotated


def get_polygon_bounds(polygon):
    # 最小包络矩阵
    if polygon is None or len(polygon) < 3:
        return None

    xmax = polygon[0]['x']
    xmin = polygon[0]['x']
    ymax = polygon[0]['y']
    ymin = polygon[0]['y']

    for point in polygon:
        if point['x'] > xmax:
            xmax = point['x']
        elif point['x'] < xmin:
            xmin = point['x']
        if point['y'] > ymax:
            ymax = point['y']
        elif point['y'] < ymin:
            ymin = point['y']

    return {
        'x': xmin,
        'y': ymin,
        'width': xmax - xmin,
        'height': ymax - ymin
    }
