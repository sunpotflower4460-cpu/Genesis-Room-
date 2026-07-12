#!/usr/bin/env python3
"""渦-依頼C 補助: 渦線の集合から「関係の複体」を組み、その次元・2-cell数を測る。

過去に確立された発見「空間の次元は"辺(関係)"でなく"第三(2-cell=三角形)"に宿る」に照らす:
渦線(1-cell)の集合に、動力学から創発した「関係(bond)」を辺として張り、3つが相互に関係すると
三角形(2-cell)が閉じる。この単体複体の次元(グラフ測地球の成長 N(r)~r^d)と2-cell数を測る。

規律: 三角形/橋は初期条件に置かない。渦の位置と符号だけ入れ、bond(関係)は「芯どうしが束縛
して近接を保つ」動力学の出力として定義し、2-cell(相互に関係する三つ組)が閉じるか・複体の
次元が1→2へ立つかを測る。決定的対照はdipole gas(辺だけ、2-cellゼロ、次元1)。
"""

import numpy as np


def bond_edges(points, bond_len):
    """points: [(x,y), ...]。距離 < bond_len の対を bond(辺) とする。"""
    n = len(points)
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            d = np.hypot(points[i][0] - points[j][0], points[i][1] - points[j][1])
            if d < bond_len:
                edges.append((i, j))
    return edges


def triangles_from_edges(n, edges):
    """辺集合から、3頂点が相互に辺で結ばれた三つ組(2-cell)を数える。"""
    adj = {i: set() for i in range(n)}
    for (i, j) in edges:
        adj[i].add(j)
        adj[j].add(i)
    tris = []
    for (i, j) in edges:
        for k in adj[i] & adj[j]:
            if k > j:  # i<j<k で重複回避（i<jは辺の順序で保証されない場合があるので下で保証）
                tri = tuple(sorted((i, j, k)))
                tris.append(tri)
    return sorted(set(tris))


def graph_ball_dimension(n, edges, max_hops=6):
    """1-skeleton上のグラフ測地球の成長 N(r) から次元を推定。
    各ノードを中心にBFSでホップ距離r以内のノード数を数え、平均 <N(r)> の
    log-log 傾きを次元とする。三角メッシュ~2, 鎖/マッチング~1。"""
    if n == 0:
        return {"ball_dim": None, "growth": []}
    adj = {i: set() for i in range(n)}
    for (i, j) in edges:
        adj[i].add(j)
        adj[j].add(i)
    # 各中心からBFS
    max_r = max_hops
    counts = np.zeros(max_r + 1)  # 累積 <N(<=r)>
    for src in range(n):
        dist = {src: 0}
        frontier = [src]
        d = 0
        while frontier and d < max_r:
            d += 1
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = d
                        nxt.append(v)
            frontier = nxt
        for r in range(max_r + 1):
            counts[r] += sum(1 for v, dd in dist.items() if dd <= r)
    counts /= n  # 平均 <N(<=r)>
    # log-log 傾き: r>=1 で N が飽和する前の範囲を使う
    rs, ns = [], []
    for r in range(1, max_r + 1):
        if counts[r] > counts[r - 1] + 1e-9 and counts[r] < 0.9 * n:
            rs.append(r)
            ns.append(counts[r])
    if len(rs) >= 2:
        slope = float(np.polyfit(np.log(rs), np.log(ns), 1)[0])
    else:
        slope = None
    return {"ball_dim": round(slope, 3) if slope is not None else None,
            "growth": [(r, round(float(counts[r]), 3)) for r in range(max_r + 1)]}


def correlation_dimension(points):
    """点集合の相関次元 C(r)~r^d (Grassberger-Procaccia)。ペア距離の累積分布の
    log-log 傾き。2Dの塊~2, 直線~1。埋め込み(位置)の幾何次元(honest floor)を測る補助指標。"""
    n = len(points)
    if n < 4:
        return None
    pts = np.array(points, dtype=float)
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            dists.append(np.hypot(pts[i, 0] - pts[j, 0], pts[i, 1] - pts[j, 1]))
    dists = np.sort(np.array(dists))
    if dists[-1] <= 0:
        return None
    # r を中間レンジ(10%..60%ile)でとり C(r)=#pairs<r/total の傾き
    lo, hi = np.percentile(dists, 10), np.percentile(dists, 60)
    rs = np.linspace(lo, hi, 8)
    rs = rs[rs > 0]
    cr = np.array([np.mean(dists < r) for r in rs])
    good = cr > 0
    if good.sum() < 3:
        return None
    slope = float(np.polyfit(np.log(rs[good]), np.log(cr[good]), 1)[0])
    return round(slope, 3)


def complex_summary(points, bond_len, max_hops=6):
    """点集合+bond距離から複体を組み、頂点/辺/2-cell/オイラー標数/球次元を返す。"""
    n = len(points)
    edges = bond_edges(points, bond_len)
    tris = triangles_from_edges(n, edges)
    euler = n - len(edges) + len(tris)  # V - E + F
    ball = graph_ball_dimension(n, edges, max_hops=max_hops)
    # 連結成分数
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for (i, j) in edges:
        parent[find(i)] = find(j)
    n_comp = len(set(find(i) for i in range(n))) if n else 0
    return {"n_vertices": n, "n_edges": len(edges), "n_2cells": len(tris),
            "triangles_per_vertex": round(len(tris) / n, 3) if n else 0.0,
            "euler_char": euler, "n_components": n_comp,
            "ball_dim": ball["ball_dim"], "ball_growth": ball["growth"],
            "correlation_dim": correlation_dimension(points)}


if __name__ == "__main__":
    # セルフテスト: 三角格子パッチ ~ dim2, 直線鎖 ~ dim1, dipole gas(離れた辺) ~ 2-cellゼロ
    def tri_lattice(nx, ny, a=1.0):
        pts = []
        for iy in range(ny):
            for ix in range(nx):
                pts.append((ix * a + (iy % 2) * a * 0.5, iy * a * np.sqrt(3) / 2))
        return pts
    print("triangular lattice 6x6:", complex_summary(tri_lattice(6, 6), bond_len=1.2))
    chain = [(i * 1.0, 0.0) for i in range(20)]
    print("1D chain:", complex_summary(chain, bond_len=1.2))
    dipoles = []
    for k in range(8):
        dipoles += [(k * 5.0, 0.0), (k * 5.0 + 1.0, 0.0)]  # 8個の離れたペア
    print("dipole gas:", complex_summary(dipoles, bond_len=1.2))
