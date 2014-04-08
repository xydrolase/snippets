#!/usr/bin/env python

__doc__ = """
Convert a correlation matrix (see example below) into a DOT file which can be 
visualized by GraphViz (more specifically, the `circo` algorithm).
"""

"""
#              v1          v2          v3         hv4        v3ab         lo1         lo2         loc        hmtp 
v1          0.0         1.5000     -0.4743      0.3128     -0.2825      0.0         0.0         0.0        -0.1746
v2         -0.2929      0.0         1.3891     -0.1970      0.0         0.0         0.0         0.0         0.0   
v3          0.0         0.0         0.0         0.8524     -0.2122      0.4098      0.0        -0.1567      0.0   
hv4         0.0         0.0         0.0         0.0        -0.1771      0.0         0.3243      0.0         0.6512
v3ab        0.0         0.0         1.3317      0.0         0.0        -0.4871      0.0         0.0         0.0   
lo1         0.0         0.0         0.0         0.0         0.0         0.0         0.5349      0.2373      0.0   
lo2         0.0         0.0         0.0         0.0         0.0         0.0         0.0         0.4579      0.0   
loc         0.0         0.0         0.0         0.0         0.6640      0.0         0.0         0.0         0.0   
hmtp        0.0         0.0         0.0         0.0         0.0         0.6069      0.2987      0.0         0.0   
"""

from argparse import ArgumentParser

import re

def color_map_discrete(val, cstart, cend, vmin, vmax):
    if val > 0.0:
        return ((hex(int(c*255))[2:]).rjust(2, '0') for c in cstart)
    else:
        return ((hex(int(c*255))[2:]).rjust(2, '0') for c in cend)

def color_map_continuous(val, cstart, cend, vmin, vmax):
    norm_val = (val - vmin) / (vmax - vmin)
    return (
        (hex(int(abs(l + norm_val * (u-l)) * 255))[2:]).rjust(2, '0')
        for l, u in zip(cstart, cend)
    )

def thickness_map(val, vmax, tl, tu):
    return tl + (tu - tl) * abs(val) / vmax

def read_matrix(fname):
    with open(fname) as f:
        colnames = re.split(r' +', f.readline().strip())[1:]
        rows = [re.split(r' +', line.strip()) for line in f if line.strip()]

    rownames = [r[0] for r in rows]
    matrix = [map(float, r[1:]) for r in rows]
    
    return rownames, colnames, matrix

def edge_to_dot(edge, padding=""):
    node_from, node_to, color, t, val, lbl = edge
    invis = ",style=invis" if val == 0.0 else ""
    if lbl is None:
        return '{4}{0} -> {1} [color="#{2}",penwidth={3}{5}]'.format(
            node_from, node_to, ''.join(map(str, color)), t, padding, invis)
    else:
        return '{5}{0} -> {1} [color="#{2}",penwidth={3},label="{4}"{6}]'.format(
            node_from, node_to, ''.join(map(str, color)), t, lbl, 
            padding, invis)

def main():
    parser = ArgumentParser(
        description="Convert a matrix into DOT language.")
    parser.add_argument("-L", "--no-label", dest="no_label",
                        default=False, action="store_true",
                        help="Suppress edge labeling.")
    parser.add_argument("-C", dest="continuous_color",
                        default=False, action="store_true",
                        help="Use continuous color mapping.")
    parser.add_argument("-T", dest="thickness", nargs=2,
                        default=(0.1, 5),
                        help="Range for line thickness.")
    parser.add_argument("-n", dest="name", default="M",
                        help="Matrix name in DOT output.")
    parser.add_argument("matrix", 
                        help="A text file containing matrix information.")
    args = parser.parse_args()

    rnames, cnames, matrix = read_matrix(args.matrix)

    maxval = max([max(r) for r in matrix])
    minval = min([min(r) for r in matrix])
    maxabs = max([abs(maxval), abs(minval)])

    cmap = color_map_continuous \
            if args.continuous_color else color_map_discrete

    edges = []
    for r in range(len(matrix)):
        m_slice = matrix[r]
        for c in range(len(m_slice)):
            if True or m_slice[c] != 0.0:
                edges.append(
                    edge_to_dot(
                        (rnames[r], cnames[c],
                         cmap(
                             m_slice[c], 
                             (1.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                             -maxabs, maxabs),
                         thickness_map(m_slice[c], maxabs, 
                                       args.thickness[0], args.thickness[1]),
                         m_slice[c],
                        str(m_slice[c]) if not args.no_label else None
                        ),
                        "    "
                    )
                )

    print """digraph %s {\n%s\n}""" % (args.name, "\n".join(edges))

if __name__ == "__main__":
    main()
