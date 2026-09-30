"""环状 DNA 条码拼接核心算法。

把每条长度 L = k+1 的读数看成德布鲁因多重图中的一条有向边
    s[:-1]  ->  s[1:]
（序列每次出现都是一条独立平行边，即一份证据）。
"每条序列恰好使用一次、相邻重叠 k = L-1 的闭环拼接" 等价于该多重图的
欧拉回路。

等价归一化（用于判断是否为同一个可信条码）：
  * 环的循环移位；
  * 整条环的反向互补；
  * 完全相同的重复读数（平行边）互换不算新答案。

枚举在固定起点上进行带剪枝的欧拉回路 DFS，同一顶点上标签相同的候选边
只取当前未用的最小输入序号一条（稳定分配），每一步都丢掉"跨过桥"的走法，
找到 2 个不同规范等价类即可判定为歧义。
"""
from __future__ import annotations

from dataclasses import dataclass

ALPHABET = frozenset("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")

MIN_READS, MAX_READS = 6, 24
MIN_LEN, MAX_LEN = 3, 8


class ValidationError(ValueError):
    """输入不满足 6~24 条、等长 3~8、仅含 ACGT 等约束。"""


@dataclass(frozen=True)
class Edge:
    u: str       # 起点 k-mer
    v: str       # 终点 k-mer
    label: str   # 原始读数
    idx: int     # 从 0 开始的输入序号


def reverse_complement(s: str) -> str:
    return s.translate(_COMP)[::-1]


def canonical_circle(seq: str) -> str:
    """环序列在"循环移位 + 整条反向互补"下的规范代表（字典序最小者）。"""
    n = len(seq)
    best = seq
    s = seq
    for _ in range(n - 1):
        s = s[1:] + s[0]
        if s < best:
            best = s
    s = reverse_complement(seq)
    for _ in range(n):
        if s < best:
            best = s
        s = s[1:] + s[0]
    return best


def validate(sequences: list[str]) -> int:
    """校验输入，返回 k = L-1；不合法抛 ValidationError。"""
    if not isinstance(sequences, list) or not sequences:
        raise ValidationError("sequences 必须是非空数组")
    if not (MIN_READS <= len(sequences) <= MAX_READS):
        raise ValidationError(
            f"序列条数必须在 {MIN_READS}~{MAX_READS} 之间，实际 {len(sequences)} 条"
        )
    norm: list[str] = []
    for i, s in enumerate(sequences):
        if not isinstance(s, str):
            raise ValidationError(f"第 {i + 1} 条序列不是字符串")
        t = s.strip().upper()
        if not (MIN_LEN <= len(t) <= MAX_LEN):
            raise ValidationError(
                f"第 {i + 1} 条序列长度必须在 {MIN_LEN}~{MAX_LEN}，实际 {len(t)}"
            )
        bad = sorted(set(t) - ALPHABET)
        if bad:
            raise ValidationError(f"第 {i + 1} 条序列含非法碱基: {''.join(bad)}")
        norm.append(t)
    lengths = {len(s) for s in norm}
    if len(lengths) != 1:
        raise ValidationError(f"所有序列必须等长，实际长度集合 {sorted(lengths)}")
    # 规范化后写回原列表（大写、去空白）
    sequences[:] = norm
    return len(norm[0]) - 1


def _components(edges: list[Edge]) -> list[list[str]]:
    """非零度顶点的弱连通分量（并查集）。"""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for e in edges:
        union(e.u, e.v)
    groups: dict[str, list[str]] = {}
    for v in sorted(parent):
        groups.setdefault(find(v), []).append(v)
    return [sorted(g) for g in groups.values()]


def _remaining_connected(mask: int, edges: list[Edge], cur: str, start: str) -> bool:
    """剩余边的所有端点（连同 cur、start）是否处于同一个弱连通分量。

    不满足说明当前走法跨过了桥，余下的边不可能再被一条连续轨迹走完并回到
    起点，可安全剪枝。
    """
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    m = mask
    while m:
        b = m & -m
        i = b.bit_length() - 1
        e = edges[i]
        ra, rb = find(e.u), find(e.v)
        if ra != rb:
            parent[rb] = ra
        m ^= b

    root = find(cur)
    if find(start) != root:
        return False
    for v in parent:
        if find(v) != root:
            return False
    return True


def _barcode_from_path(path: list[int], edges: list[Edge], k: int) -> str:
    """由边序列拼出长度为 n 的环状条码（闭环后末尾 k-mer 与起点重合）。"""
    s = edges[path[0]].label
    for eid in path[1:]:
        s += edges[eid].label[-1]
    return s[: len(path)]  # 丢掉与起点重合的末尾 k 个字符


def _align_order(
    path: list[int], circle: str, canon: str
) -> tuple[list[int], str]:
    """尽量用纯循环移位把使用次序对齐到规范条码；否则只能经反向互补归一。"""
    n = len(circle)
    s = circle
    for r in range(n):
        if s == canon:
            return path[r:] + path[:r], "cyclic_rotation" if r else "identical"
        s = s[1:] + s[0]
    return path, "reverse_complement"


def _enumerate_classes(
    edges: list[Edge], start: str, k: int, limit: int = 2
) -> list[tuple[str, list[int], str]]:
    """枚举至多 limit 个不同规范等价类，返回 (规范条码, 边序号路径, 归一关系)。"""
    m = len(edges)
    full = (1 << m) - 1

    outgoing: dict[str, list[int]] = {}
    for e in edges:
        outgoing.setdefault(e.u, []).append(e.idx)

    found: list[tuple[str, list[int], str]] = []

    def dfs(cur: str, used: int, path: list[int]) -> None:
        if len(found) >= limit:
            return
        if used == full:
            if cur != start:
                return
            circle = _barcode_from_path(path, edges, k)
            canon = canonical_circle(circle)
            if not any(c == canon for c, _, _ in found):
                order, rel = _align_order(path, circle, canon)
                found.append((canon, order, rel))
            return

        # 同标签平行边只尝试当前未用的最小序号 -> 重复读数互换不产生歧义，
        # 且使用次序按输入序号稳定分配。
        choices: dict[str, int] = {}
        for eid in outgoing.get(cur, ()):
            if not (used >> eid) & 1:
                choices.setdefault(edges[eid].label, eid)

        can_finish = (full ^ used).bit_count() == 1
        for eid in choices.values():
            e = edges[eid]
            nused = used | (1 << eid)
            if not can_finish and not _remaining_connected(
                full ^ nused, edges, e.v, start
            ):
                continue
            dfs(e.v, nused, path + [eid])
            if len(found) >= limit:
                return

    dfs(start, 0, [])
    return found


DIRECTION_FORWARD = "forward"
DIRECTION_REVERSE = "reverse_complement"


def _evidence(
    order: list[int],
    edges: list[Edge],
    k: int,
    directions: list[str] | None = None,
) -> list[dict]:
    """相邻（含首尾相接）重叠证据。

    方向未知模式下 directions 非空，逐项补充该读数实际采用的方向与
    定向后序列（原序列或其反向互补）。
    """
    out = []
    m = len(order)
    for pos, eid in enumerate(order):
        a = edges[eid]
        b = edges[order[(pos + 1) % m]]
        overlap = a.label[1:]
        assert overlap == b.label[:-1]
        item = {
            "position": pos,
            "prev": a.idx + 1,
            "next": b.idx + 1,
            "overlap": overlap,
            "overlap_length": k,
            "appended_base": b.label[-1],
        }
        if directions is not None:
            item["orientation"] = directions[a.idx]
            item["oriented_sequence"] = a.label
        out.append(item)
    return out


def _witness(
    canon: str,
    order0: list[int],
    relation: str,
    edges: list[Edge],
    k: int,
    directions: list[str] | None = None,
) -> dict:
    return {
        "canonical_barcode": canon,
        "barcode": _barcode_from_path(order0, edges, k),
        "canonical_relation": relation,
        "order": [i + 1 for i in order0],
        "evidence": _evidence(order0, edges, k, directions),
    }


def _assemble_known(seq: list[str], k: int) -> dict:
    """方向已知（默认模式）：固定方向建图、判定并生成见证。"""
    edges = [
        Edge(s[:-1], s[1:], s, i)
        for i, s in enumerate(seq)
    ]

    indeg: dict[str, int] = {}
    outdeg: dict[str, int] = {}
    for e in edges:
        outdeg[e.u] = outdeg.get(e.u, 0) + 1
        indeg[e.v] = indeg.get(e.v, 0) + 1
    vertices = sorted(set(indeg) | set(outdeg))

    reasons: list[dict] = []

    imbalance = [
        {"vertex": v, "in_degree": indeg.get(v, 0), "out_degree": outdeg.get(v, 0)}
        for v in vertices
        if indeg.get(v, 0) != outdeg.get(v, 0)
    ]
    if imbalance:
        reasons.append(
            {
                "code": "degree_imbalance",
                "message": "存在入度不等于出度的 k-mer，无法形成闭环",
                "vertices": imbalance,
            }
        )

    comps = _components(edges)
    if len(comps) > 1:
        reasons.append(
            {
                "code": "fragmented_graph",
                "message": "非零度顶点落在多个互不连通的分量中，无法拼成单一闭环",
                "component_count": len(comps),
                "components": comps,
            }
        )

    if reasons:
        return {
            "status": "no_solution",
            "overlap_length": k,
            "read_count": len(edges),
            "reasons": reasons,
        }

    start = vertices[0]  # 回路经过所有顶点，固定最小顶点作为枚举起点
    classes = _enumerate_classes(edges, start, k, limit=2)

    if len(classes) == 1:
        canon, order0, rel = classes[0]
        return {
            "status": "unique",
            "overlap_length": k,
            "read_count": len(edges),
            **_witness(canon, order0, rel, edges, k),
            "normalizations": ["cyclic_rotation", "reverse_complement"],
        }

    if not classes:
        # 度数平衡且连通时必然存在欧拉回路；到不了这里只是防御性兜底。
        return {
            "status": "no_solution",
            "overlap_length": k,
            "read_count": len(edges),
            "reasons": [
                {
                    "code": "euler_search_failed",
                    "message": "图满足必要条件但未找到欧拉回路（内部错误）",
                }
            ],
        }

    witnesses = [_witness(c, o, r, edges, k) for c, o, r in classes]
    return {
        "status": "ambiguous",
        "overlap_length": k,
        "read_count": len(edges),
        "message": "存在多个不同规范等价类的闭环拼法，给出两条不同规范见证",
        "witnesses": witnesses,
        "normalizations": ["cyclic_rotation", "reverse_complement"],
    }


def _folded_components(seq: list[str], k: int) -> list[list[str]]:
    """把每个 k-mer 与它的反向互补折叠成超级节点后的非空弱连通分量。

    读数 s 的两个端点为 s[:-1]、s[1:]；取反向互补后端点变为
    rc(s[1:])、rc(s[:-1])，折叠（w~rc(w)）之后是同一条无向边。因此该
    折叠图与逐条定向【无关】：折叠图不连通则任何定向都不可能得到单一
    欧拉闭环（必要条件）。
    """
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for s in seq:
        u, v = s[:-1], s[1:]
        fu = min(u, reverse_complement(u))
        fv = min(v, reverse_complement(v))
        union(fu, fv)

    groups: dict[str, list[str]] = {}
    for node in sorted(parent):
        groups.setdefault(find(node), []).append(node)
    return [sorted(members) for members in groups.values()]


# ---------------------------------------------------------------------------
# 方向未知模式：每条读数【独立】可按原序列或其反向互补参与拼接，但每个输入
# 序号只能使用一次。注意 s 与 rc(s) 是两条不同的读数，定向变量也相互独立
# （环上完全可能同时出现 w 与 rc(w)，不能强制它们同向）；只有完全相同的
# 平行读数可互换，用"该标签取反向互补的份数 x ∈ 0..c"压缩搜索。
#
# 定向与闭环次序【联合】求解（不是整批翻转两次尝试）：
#   1) 在"每个标签反向份数"变量上回溯，增量维护各 k-mer 的 出度-入度 差；
#      用剩余变量在该 k-mer 上可达贡献的逐点上下界剪枝（平衡是闭环必要
#      条件）；
#   2) 完整定向再检查弱连通并枚举欧拉回路、按规范类去重，找到 2 个不同
#      规范类即停。
# 全局对称剪枝：所有读数整体翻转只把环换成其反向互补（同一规范类），故把
# 第 1 条读数的定向标签固定为 min(s, rc(s)) 不丢解，搜索减半且结果稳定。
# ---------------------------------------------------------------------------


def _assemble_unknown(seq: list[str], k: int) -> dict:
    n = len(seq)

    # 定向无关的必要条件预检：k-mer 反向互补折叠图必须连通。
    folded = _folded_components(seq, k)
    if len(folded) > 1:
        return {
            "status": "no_solution",
            "orientation_mode": "unknown",
            "overlap_length": k,
            "read_count": n,
            "reasons": [
                {
                    "code": "orientation_fragmented_graph",
                    "message": (
                        "方向未知：把每个 k-mer 与其反向互补折叠为同一节点后，"
                        "读数端点仍落在多个互不连通的分量中。读数无论正向还是反向"
                        "互补，在折叠图上都是同一条边，故该断裂与逐条定向无关，"
                        "任何定向都无法拼成单一闭环（不能用固定方向下的度数失衡"
                        "解释）。"
                    ),
                    "component_count": len(folded),
                    "folded_components": folded,
                }
            ],
        }

    # 按原始标签分组（平行读数）；组序按首次出现，保证稳定
    groups: list[dict] = []
    group_of: dict[str, int] = {}
    for i, s in enumerate(seq):
        if s in group_of:
            groups[group_of[s]]["indices"].append(i)
        else:
            group_of[s] = len(groups)
            rc = reverse_complement(s)
            groups.append(
                {
                    "label": s,
                    "rc": rc,
                    "indices": [i],
                    "fixed": False,
                }
            )
    gcount_all = len(groups)
    for g in groups:
        g["count"] = len(g["indices"])
        g["palindrome"] = g["label"] == g["rc"]

    # 回文读数（s==rc(s)）定向不改变图，不是真正的变量，直接固定为正向。
    variable = [g for g in range(gcount_all) if not groups[g]["palindrome"]]
    if variable:
        # 整批翻转对称：固定一个非回文组中代表读数的定向标签为 min(s, rc)
        fg = group_of[seq[0]]
        fixed_group = fg if not groups[fg]["palindrome"] else variable[0]
        groups[fixed_group]["fixed"] = True
    else:
        fixed_group = -1
    gcount = len(variable)

    # 每组取 x 份反向互补（x∈[0,c]）时，对各 k-mer 的 出度-入度 贡献：
    #   x=0：+c@u0, -c@v0；每多一份反向：u0 -1、v0 +1、u1 +1、v1 -1。
    # 偶数 k 时 k-mer 自身可能是回文，导致 u0==v1 / v0==u1，故四个角色
    # 对同一节点的系数必须【相加】而不是覆盖。
    for g in groups:
        s, rc, c = g["label"], g["rc"], g["count"]
        u0, v0 = s[:-1], s[1:]
        u1, v1 = rc[:-1], rc[1:]
        coef: dict[str, int] = {}
        for node, wgt in ((u0, c), (v0, -c)):
            coef[node] = coef.get(node, 0) + wgt
        step: dict[str, int] = {}
        for node, wgt in ((u0, -1), (v0, 1), (u1, 1), (v1, -1)):
            step[node] = step.get(node, 0) + wgt
        g["coef"] = coef
        g["step"] = step
        # x ∈ [0, c] 时各节点贡献的最小值/最大值
        bounds: dict[str, tuple[int, int]] = {}
        for node in set(coef) | set(step):
            vals = [
                coef.get(node, 0) + x * step.get(node, 0)
                for x in range(c + 1)
            ]
            bounds[node] = (min(vals), max(vals))
        g["bounds"] = bounds

    # 固定组先赋值（对称性），其余变量组按"份数多、标签大"优先使剪枝尽早生效
    ordered = [fixed_group] + sorted(
        (g for g in variable if g != fixed_group),
        key=lambda g: (-groups[g]["count"], groups[g]["label"]),
    )

    all_nodes = sorted({node for g in groups for node in g["bounds"]})

    # 后缀（剩余组）对各节点贡献的上下界之和
    suffix_lo: list[dict[str, int]] = [dict() for _ in range(gcount + 1)]
    suffix_hi: list[dict[str, int]] = [dict() for _ in range(gcount + 1)]
    for pos in range(gcount - 1, -1, -1):
        g = groups[ordered[pos]]
        for node in all_nodes:
            lo, hi = g["bounds"].get(node, (0, 0))
            suffix_lo[pos][node] = suffix_lo[pos + 1].get(node, 0) + lo
            suffix_hi[pos][node] = suffix_hi[pos + 1].get(node, 0) + hi

    found: list[tuple[str, tuple[int, ...], list[int], str]] = []
    found_canons: set[str] = set()
    stats = {"assignments_tested": 0, "disconnected": 0}
    # 有效搜索空间：固定组（非回文）的锚定方向钉死，取值数为 c；其余组 c+1
    search_space = 1
    for pos in range(gcount):
        gi = ordered[pos]
        search_space *= groups[gi]["count"] if groups[gi]["fixed"] else groups[gi]["count"] + 1

    def build(xs: tuple[int, ...]) -> tuple[list[Edge], list[str]]:
        """xs 为按 ordered 排列的各组"反向副本总数"；构图与逐条方向。

        普通组：反向互补稳定地分给输入序号最小的 x 份。固定组：锚定副本
        （最小序号）方向钉死为 min(label, rc)；其余 x 或 x-1 份反向分给
        非锚定副本中序号最大者，从而每个全局翻转轨道恰好出现一次。
        """
        rev_count = {ordered[pos]: x for pos, x in enumerate(xs)}
        dirs = [DIRECTION_FORWARD] * n
        edges: list[Edge] = [None] * n  # type: ignore[list-item]
        for gi, g in enumerate(groups):
            x = rev_count.get(gi, 0)
            indices = g["indices"]
            anchor, others = indices[0], indices[1:]
            if g["fixed"]:
                anchor_reversed = reverse_complement(g["label"]) < g["label"]
                extra = x - (1 if anchor_reversed else 0)
                rev_set = set(others[len(others) - extra:]) if extra else set()
                if anchor_reversed:
                    rev_set.add(anchor)
            else:
                rev_set = set(indices[:x])
            for i in indices:
                if i in rev_set:
                    s2 = reverse_complement(seq[i])
                    dirs[i] = DIRECTION_REVERSE
                else:
                    s2 = seq[i]
                    dirs[i] = DIRECTION_FORWARD
                edges[i] = Edge(s2[:-1], s2[1:], s2, i)
        return edges, dirs

    def evaluate(xst: tuple[int, ...]) -> None:
        """完整定向：平衡（调用方已保证）→ 连通 → 欧拉枚举 → 规范类去重。"""
        stats["assignments_tested"] += 1
        edges, _dirs = build(xst)
        if len(_components(edges)) > 1:
            stats["disconnected"] += 1
            return
        vertices = sorted({e.u for e in edges} | {e.v for e in edges})
        limit = 2 - len(found)
        classes = _enumerate_classes(edges, vertices[0], k, limit=limit)
        for canon, order0, rel in classes:
            if canon in found_canons:
                continue
            found_canons.add(canon)
            found.append((canon, xst, order0, rel))
            if len(found) >= 2:
                return

    def backtrack(pos: int, xs: list[int], diff: dict[str, int]) -> None:
        if len(found) >= 2:
            return
        g = groups[ordered[pos]]
        # 固定组：锚定副本（最小序号）方向钉死为 min(label, rc)，变量 x 是
        # 反向副本【总数】。锚定正向时 x ∈ 0..c-1，锚定反向时 x ∈ 1..c；
        # 每个全局翻转对（r 与 c-r）恰好出现一次。
        if g["fixed"]:
            c0 = g["count"]
            anchor = g["label"]
            if reverse_complement(anchor) < anchor:
                x_values = range(1, c0 + 1)   # 锚定反向
            else:
                x_values = range(0, c0)       # 锚定正向
        else:
            x_values = range(0, g["count"] + 1)

        touched = set(g["coef"]) | set(g["step"])
        for x in x_values:
            for node in touched:
                diff[node] = (
                    diff.get(node, 0)
                    + g["coef"].get(node, 0)
                    + x * g["step"].get(node, 0)
                )

            nxt = pos + 1
            # 终值需为 0：剩余组在该节点的贡献 ∈ [lo,hi]，故 -dv 必须落在此区间
            feasible = all(
                -suffix_hi[nxt].get(node, 0) <= dv <= -suffix_lo[nxt].get(node, 0)
                for node, dv in diff.items()
            )
            if feasible and nxt < gcount:
                xs.append(x)
                backtrack(nxt, xs, diff)
                xs.pop()
            elif feasible:
                # 末层边界 lo=hi=0，feasible 已保证所有节点 diff==0
                evaluate(tuple(list(xs) + [x]))

            for node in touched:
                diff[node] -= (
                    g["coef"].get(node, 0) + x * g["step"].get(node, 0)
                )

    if gcount:
        # 回文组不是变量，但其定向（固定正向）仍对度数平衡有恒定贡献，
        # 作为回溯基线计入 diff；它们不出现在任何 suffix 边界中。
        baseline: dict[str, int] = {}
        for gi, g in enumerate(groups):
            if g["palindrome"]:
                for node, val in g["coef"].items():
                    baseline[node] = baseline.get(node, 0) + val
        backtrack(0, [], baseline)
    else:
        # 全部读数均为回文（s==rc(s)）：定向不影响图，直接评估一次
        evaluate(())

    if not found:
        group_details = [
            {
                "label": g["label"],
                "reverse_complement": g["rc"],
                "multiplicity": g["count"],
                "reversed_copies_domain": [0, g["count"]],
            }
            for g in groups
        ]
        return {
            "status": "no_solution",
            "orientation_mode": "unknown",
            "overlap_length": k,
            "read_count": n,
            "reasons": [
                {
                    "code": "orientation_constraints_unsatisfiable",
                    "message": (
                        "方向未知：对每条读数的正向/反向互补定向与闭环次序进行了"
                        "联合穷举（相同平行读数只按反向份数计；所有读数整体翻转只"
                        "产生同一规范类，故把首个非回文标签的锚定副本定向固定为其"
                        "两方向的字典序较小者），不存在任何逐条定向能同时满足各 "
                        "k-mer 出入度平衡、弱连通并形成使用全部读数恰好一次的欧拉"
                        "闭环。注意：固定方向下的度数失衡在此模式下不成立、不能"
                        "作为结论。"
                    ),
                    "orientation_groups": group_details,
                    "group_assignments_evaluated": stats["assignments_tested"],
                    "disconnected_balanced_assignments": stats["disconnected"],
                    "raw_group_assignment_space": search_space,
                    "necessary_conditions": [
                        "每个 k-mer 的入度等于出度",
                        "全部非零度顶点处于同一弱连通分量",
                        "每个输入序号恰好使用一次",
                    ],
                    "review_hint": (
                        "可按 orientation_groups 复核：对每组独立取 0..multiplicity "
                        "份反向互补，全部组合均无法通过上述必要条件并欧拉闭环"
                    ),
                }
            ],
        }

    def witness_unknown(
        canon: str, xst: tuple[int, ...], order0: list[int], rel: str
    ) -> dict:
        edges, dirs = build(xst)
        w = _witness(canon, order0, rel, edges, k, dirs)
        w["read_orientations"] = [
            {
                "read": i + 1,
                "original_sequence": seq[i],
                "orientation": dirs[i],
                "oriented_sequence": (
                    seq[i] if dirs[i] == DIRECTION_FORWARD
                    else reverse_complement(seq[i])
                ),
            }
            for i in range(n)
        ]
        return w

    common = {
        "orientation_mode": "unknown",
        "overlap_length": k,
        "read_count": n,
        "normalizations": ["cyclic_rotation", "reverse_complement"],
        "orientation_note": (
            "整批反向互补只产生同一规范条码类别，多种等价定向不构成歧义"
        ),
    }

    if len(found) == 1:
        canon, xst, order0, rel = found[0]
        return {
            "status": "unique",
            **common,
            **witness_unknown(canon, xst, order0, rel),
        }

    return {
        "status": "ambiguous",
        **common,
        "message": "存在多个不同规范等价类的定向闭环拼法，给出两条不同规范见证",
        "witnesses": [witness_unknown(*f) for f in found],
    }


def assemble(sequences: list[str], unknown_orientation: bool = False) -> dict:
    """主入口：返回唯一 / 歧义 / 无解三类结果之一。

    unknown_orientation=False（默认）时保持既有输入、结论与证据结构不变；
    为 True 时启用方向未知模式，逐条读数的定向与闭环次序联合求解。
    """
    seq = list(sequences)
    k = validate(seq)
    if unknown_orientation:
        return _assemble_unknown(seq, k)
    return _assemble_known(seq, k)