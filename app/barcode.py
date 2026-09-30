"""环状 DNA 条码拼接核心算法。

把每条长度 L = k+1 的读数看成德布鲁因多重图中的一条有向边
    s[:-1]  ->  s[1:]
（序列每次出现都是一条独立平行边，即一份证据）。
"每条序列恰好使用一次、相邻重叠 k = L-1 的闭环拼接" 等价于该多重图的
欧拉回路。

两种模式：

* 固定方向（默认）：每条读数只能按提交时的方向参与拼接。
* 方向未知（unknown_orientation=True）：部分测序批次丢失了单条读数的链
  方向，每条读数可独立按原序列或其反向互补参与拼接，但每个输入序号仍只能
  使用一次；"逐条定向"与"闭环次序"在同一次 DFS 中联合求解（不能退化为
  只尝试整批翻转）。

等价归一化（用于判断是否为同一个可信条码）：
  * 环的循环移位；
  * 整条环的反向互补（方向未知模式下每条读数的采用方向随之翻转）；
  * 完全相同的重复读数（平行边）互换不算新答案。

枚举在固定起点上进行带剪枝的欧拉回路 DFS，同一顶点上定向后标签相同的
候选边只取当前未用的最小输入序号一条（稳定分配），每一步都丢掉"跨过桥"
的走法，找到 2 个不同规范等价类即可判定为歧义。
"""
from __future__ import annotations

from dataclasses import dataclass

ALPHABET = frozenset("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")

MIN_READS, MAX_READS = 6, 24
MIN_LEN, MAX_LEN = 3, 8

# 方向未知模式穷尽枚举时允许访问的 DFS 状态上限（远超约束下的真实需求，
# 仅为防御异常输入导致的指数膨胀）。
SEARCH_STATE_LIMIT = 2_000_000


class ValidationError(ValueError):
    """输入不满足 6~24 条、等长 3~8、仅含 ACGT 等约束。"""


@dataclass(frozen=True)
class Edge:
    u: str       # 起点 k-mer
    v: str       # 终点 k-mer
    label: str   # 原始读数
    idx: int     # 从 0 开始的输入序号


@dataclass(frozen=True)
class OrientedEdge:
    u: str           # 起点 k-mer
    v: str           # 终点 k-mer
    label: str       # 定向后的读数
    orientation: str # "+" 原序列 / "-" 反向互补


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


# ---------------------------------------------------------------------------
# 固定方向模式
# ---------------------------------------------------------------------------

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


def _barcode_from_labels(labels: list[str]) -> str:
    """由定向后读数序列拼出长度为 n 的环状条码（闭环后末尾 k-mer 与起点重合）。"""
    s = labels[0]
    for label in labels[1:]:
        s += label[-1]
    return s[: len(labels)]  # 丢掉与起点重合的末尾 k 个字符


def _barcode_from_path(path: list[int], edges: list[Edge], k: int) -> str:
    """固定方向模式的条码拼接包装。"""
    return _barcode_from_labels([edges[eid].label for eid in path])


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


def _evidence(order: list[int], edges: list[Edge], k: int) -> list[dict]:
    """相邻（含首尾相接）重叠证据。"""
    out = []
    m = len(order)
    for pos, eid in enumerate(order):
        a = edges[eid]
        b = edges[order[(pos + 1) % m]]
        overlap = a.label[1:]
        assert overlap == b.label[:-1]
        out.append(
            {
                "position": pos,
                "prev": a.idx + 1,
                "next": b.idx + 1,
                "overlap": overlap,
                "overlap_length": k,
                "appended_base": b.label[-1],
            }
        )
    return out


# ---------------------------------------------------------------------------
# 方向未知模式
# ---------------------------------------------------------------------------

def _comp_group_key(x: str) -> tuple[str, ...]:
    """k-mer 与其反向互补构成的"配对组"规范键（回文 k-mer 自成单元素组）。"""
    y = reverse_complement(x)
    return (x,) if x == y else ((x, y) if x < y else (y, x))


def _group_precheck(seq: list[str]) -> list[dict]:
    """方向未知问题的组级必要条件（完备枚举前的快速、可复核否决）。

    一条读数正向连接 k-mer 组 (G(s[:-1]), G(s[1:]))；取反向互补时两端点
    同时换成各自的互补 k-mer，所属组不变、仅组间方向翻转。因此任何可行的
    逐条定向都在"互补配对组"的无向多重图上诱导一条欧拉回路，必要条件为：
    每个非孤立组关联边数为偶数，且非孤立组全部连通。
    不满足时，无论怎样逐条定向都不可能闭环——这比沿用固定方向下的度数
    失衡更准确（翻转单条读数会同时改变两个 k-mer 的入/出度）。
    """
    parent: dict[tuple[str, ...], tuple[str, ...]] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    degree: dict[tuple[str, ...], int] = {}
    group_reads: dict[tuple[str, ...], list[int]] = {}

    def touch(g, i: int) -> None:
        degree[g] = degree.get(g, 0) + 1
        group_reads.setdefault(g, []).append(i + 1)

    for i, s in enumerate(seq):
        ga = _comp_group_key(s[:-1])
        gb = _comp_group_key(s[1:])
        # 无向边：自环在度数中计 2
        touch(ga, i)
        touch(gb, i)
        union(ga, gb)

    reasons: list[dict] = []

    odd = [
        {"group": list(g), "degree": degree[g], "reads": group_reads[g]}
        for g in sorted(degree)
        if degree[g] % 2 == 1
    ]
    if odd:
        reasons.append(
            {
                "code": "orientation_group_imbalance",
                "message": (
                    "把 k-mer 与其反向互补归为同一配对组后，存在关联读数条数为"
                    "奇数的组；一条读数无论正向还是反向互补，其两端所属组不变，"
                    "故奇度组不可能靠逐条定向消除，闭环方向约束不可满足"
                ),
                "odd_degree_groups": odd,
            }
        )

    comps: dict[tuple[str, ...], list[tuple[str, ...]]] = {}
    for g in sorted(parent):
        comps.setdefault(find(g), []).append(g)
    components = [
        sorted((list(g) for g in groups), key=lambda m: m[0])
        for groups in comps.values()
    ]
    components.sort(key=lambda c: c[0][0])
    if len(components) > 1:
        reasons.append(
            {
                "code": "orientation_group_fragmented",
                "message": (
                    "互补配对组落在多个互不连通的分量中；逐条定向只能在组间"
                    "翻转方向，不能把不同分量连成同一个闭环"
                ),
                "component_count": len(components),
                "components": components,
            }
        )
    return reasons


def _enumerate_unknown(
    seq: list[str], k: int, limit: int = 2
) -> tuple[list[tuple[str, list[tuple[int, str]], str]], dict]:
    """联合枚举"每条读数二选一定向 + 闭环次序"。

    锚定输入序号 0 为首条边（它可正向、可反向互补）；任何闭环解都含该读数，
    循环移位不损失解。返回至多 limit 个 (规范条码, [(序号, 方向)], 归一关系)
    以及穷尽枚举的搜索统计。
    """
    n = len(seq)
    full = (1 << n) - 1

    plus: list[OrientedEdge] = []
    minus: list[OrientedEdge | None] = []
    for s in seq:
        r = reverse_complement(s)
        plus.append(OrientedEdge(s[:-1], s[1:], s, "+"))
        # 回文读数（s == rc(s)）两个方向是同一条边，只保留一份，统一记正向。
        minus.append(None if r == s else OrientedEdge(r[:-1], r[1:], r, "-"))

    # 预算：组级无向边与四端连通对（均与定向无关）
    gpair: dict[str, tuple[str, ...]] = {}
    group_edges: list[tuple[tuple[str, ...], tuple[str, ...]]] = []
    hyper_pairs: list[tuple[tuple[str, str], tuple[str, str]]] = []
    for s in seq:
        a, b = s[:-1], s[1:]
        ra, rb = reverse_complement(a), reverse_complement(b)
        ga, gb = _comp_group_key(a), _comp_group_key(b)
        for x in (a, b, ra, rb):
            gpair.setdefault(x, _comp_group_key(x))
        group_edges.append((ga, gb))
        hyper_pairs.append(((a, b), (rb, ra)))  # 正向端对 / 反向互补端对

    found: list[tuple[str, list[tuple[int, str]], str]] = []
    stats = {
        "states_visited": 0,
        "no_continuation": 0,
        "group_balance_pruned": 0,
        "endpoint_connectivity_pruned": 0,
        "end_gap": 0,
        "limit_truncated": False,
    }
    dead_ends: list[dict] = []
    seen_dead: set[tuple[str, int]] = set()

    def pick_edge(i: int, orient: str) -> OrientedEdge:
        return plus[i] if orient == "+" else minus[i]  # type: ignore[return-value]

    def choices_at(cur: str, used: int) -> list[tuple[int, str]]:
        """当前顶点可接的 (读数, 方向)，定向后同标签只留最小序号（稳定分配）。"""
        choices: dict[str, tuple[int, str]] = {}
        m = full ^ used
        while m:
            b = m & -m
            i = b.bit_length() - 1
            m ^= b
            for e in (plus[i], minus[i]):
                if e is not None and e.u == cur:
                    old = choices.get(e.label)
                    # 元组比较先比序号；同序号时 "+" < "-"，确定性优先正向
                    if old is None or (i, e.orientation) < old:
                        choices[e.label] = (i, e.orientation)
        return sorted(choices.values())

    def group_feasible(mask: int, cur: str, start: str) -> bool:
        """剩余读数（加一条 start 组到 cur 组的虚边）能否形成组级平衡迹。

        无向多重图存在连接两顶点、用完所有边的迹，当且仅当恰有 0 或 2 个
        奇度顶点且非孤立顶点连通；加虚边后统一为"全偶 + 连通"，是原问题
        的必要松弛。
        """
        parent: dict[tuple[str, ...], tuple[str, ...]] = {}

        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        degree: dict[tuple[str, ...], int] = {}

        def add(ga, gb) -> None:
            degree[ga] = degree.get(ga, 0) + 1
            degree[gb] = degree.get(gb, 0) + 1
            ra, rb = find(ga), find(gb)
            if ra != rb:
                parent[rb] = ra

        m = mask
        while m:
            b = m & -m
            i = b.bit_length() - 1
            m ^= b
            add(*group_edges[i])

        gc, gs = gpair[cur], gpair[start]
        add(gc, gs)  # 虚边（同组时为自环，度数 +2 不影响奇偶）

        if any(d % 2 == 1 for d in degree.values()):
            return False
        root = find(gc)
        return all(find(g) == root for g in degree)

    def endpoints_connected(mask: int, cur: str, start: str) -> bool:
        """端点级松弛：把每条剩余读数的 4 个候选端点并入同一个超边。

        读数 i 实际只会选用 (a,b) 或 (rb,ra) 中的一对，实际端对是该超边的
        子集；因此真实迹连通时超图必然连通，反之不连通即可安全剪枝。
        （不能只把两个端对各自并起来：未选用方向的备选端点会形成假分量。）
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

        m = mask
        while m:
            b = m & -m
            i = b.bit_length() - 1
            m ^= b
            (a1, b1), (a2, b2) = hyper_pairs[i]
            # 四个候选端点属于同一条读数（超边），全部并到一起
            for t in (b1, a2, b2):
                union(a1, t)

        root = find(cur)
        if find(start) != root:
            return False
        return all(find(v) == root for v in parent)

    def record_dead_end(cur: str, remaining: int) -> None:
        key = (cur, remaining)
        if key in seen_dead or len(dead_ends) >= 3:
            return
        seen_dead.add(key)
        candidates = []
        m = remaining
        while m:
            b = m & -m
            i = b.bit_length() - 1
            m ^= b
            for e in (plus[i], minus[i]):
                if e is not None:
                    candidates.append(
                        {
                            "read": i + 1,
                            "orientation": (
                                "forward" if e.orientation == "+"
                                else "reverse_complement"
                            ),
                            "start_kmer": e.u,
                            "oriented_sequence": e.label,
                        }
                    )
        dead_ends.append(
            {
                "at_kmer": cur,
                "remaining_reads": [
                    i + 1 for i in range(n) if (remaining >> i) & 1
                ],
                "candidate_edges": candidates,
                "note": "该 k-mer 下所有未用读数的两种定向都没有可接的边",
            }
        )

    def record(path: list[tuple[int, str]]) -> None:
        labels = [
            seq[i] if orient == "+" else reverse_complement(seq[i])
            for i, orient in path
        ]
        canon = canonical_circle(_barcode_from_labels(labels))
        if any(c == canon for c, _, _ in found):
            return
        items, rel = _normalize_unknown_path(path, seq, canon)
        found.append((canon, items, rel))

    def dfs(cur: str, start_v: str, used: int,
            path: list[tuple[int, str]]) -> None:
        stats["states_visited"] += 1
        if stats["states_visited"] > SEARCH_STATE_LIMIT:
            stats["limit_truncated"] = True
            return
        if len(found) >= limit:
            return
        if used == full:
            if cur != start_v:
                stats["end_gap"] += 1
            else:
                record(path)
            return

        remaining = full ^ used
        choices = choices_at(cur, used)
        if not choices:
            stats["no_continuation"] += 1
            record_dead_end(cur, remaining)
            return

        last = (remaining & (remaining - 1)) == 0
        for i, orient in choices:
            if len(found) >= limit or stats["limit_truncated"]:
                return
            e = pick_edge(i, orient)
            nused = used | (1 << i)
            if last:
                # 只剩一条边：组级松弛可能通过而端点对不上，精确核对闭环
                if e.v != start_v:
                    stats["end_gap"] += 1
                    continue
            else:
                rest = full ^ nused
                if not group_feasible(rest, e.v, start_v):
                    stats["group_balance_pruned"] += 1
                    continue
                if not endpoints_connected(rest, e.v, start_v):
                    stats["endpoint_connectivity_pruned"] += 1
                    continue
            dfs(e.v, start_v, nused, path + [(i, orient)])

    # 起始读数两个方向各自成锚（回文时 minus 为 None，只跑一次）
    for first in (plus[0], minus[0]):
        if first is None:
            continue
        dfs(first.v, first.u, 1, [(0, first.orientation)])
        if len(found) >= limit or stats["limit_truncated"]:
            break

    stats["dead_ends"] = dead_ends
    return found, stats


def _normalize_unknown_path(
    path: list[tuple[int, str]], seq: list[str], canon: str
) -> tuple[list[tuple[int, str]], str]:
    """把联合枚举路径对齐到规范条码。

    先尝试纯循环移位；否则把整环反向互补——新环上边次序为
    rc(e0), rc(e_{n-1}), ..., rc(e1)，每条读数的采用方向随之翻转。
    """
    n = len(path)

    def label_of(item: tuple[int, str]) -> str:
        i, orient = item
        return seq[i] if orient == "+" else reverse_complement(seq[i])

    labels = [label_of(item) for item in path]
    circle = _barcode_from_labels(labels)

    s = circle
    for r in range(n):
        if s == canon:
            return path[r:] + path[:r], (
                "identical" if r == 0 else "cyclic_rotation"
            )
        s = s[1:] + s[0]

    def flip(orient: str) -> str:
        return "-" if orient == "+" else "+"

    rev = [(path[0][0], flip(path[0][1]))]
    rev += [(idx, flip(o)) for idx, o in reversed(path[1:])]
    # 回文读数的两个方向序列相同，统一记为正向，避免见证中出现伪方向。
    rev = [
        (idx, "+") if reverse_complement(seq[idx]) == seq[idx] else (idx, o)
        for idx, o in rev
    ]
    rev_labels = [label_of(item) for item in rev]
    circle2 = _barcode_from_labels(rev_labels)
    assert canonical_circle(circle2) == canon

    s = circle2
    for r in range(n):
        if s == canon:
            return rev[r:] + rev[:r], "reverse_complement"
        s = s[1:] + s[0]
    raise AssertionError("unreachable")  # pragma: no cover


def _unknown_witness(
    canon: str, items: list[tuple[int, str]], relation: str,
    seq: list[str], k: int,
) -> dict:
    """方向未知模式的可逐项复核见证。"""
    n = len(items)
    labels = [
        seq[i] if orient == "+" else reverse_complement(seq[i])
        for i, orient in items
    ]
    order = [i + 1 for i, _ in items]
    orientations = [
        "forward" if orient == "+" else "reverse_complement"
        for _, orient in items
    ]
    evidence = []
    for pos, ((i, orient), label) in enumerate(zip(items, labels)):
        nlabel = labels[(pos + 1) % n]
        overlap = label[1:]
        assert overlap == nlabel[:-1]
        evidence.append(
            {
                "position": pos,
                "read": i + 1,
                "prev": i + 1,
                "next": items[(pos + 1) % n][0] + 1,
                "orientation": (
                    "forward" if orient == "+" else "reverse_complement"
                ),
                "original_sequence": seq[i],
                "oriented_sequence": label,
                "overlap": overlap,
                "overlap_length": k,
                "appended_base": nlabel[-1],
            }
        )
    return {
        "canonical_barcode": canon,
        "barcode": canon,  # 已按循环移位/反向互补对齐到规范代表
        "canonical_relation": relation,
        "order": order,
        "orientations": orientations,
        "evidence": evidence,
    }


def _assemble_unknown(seq: list[str], k: int) -> dict:
    base = {
        "overlap_length": k,
        "read_count": len(seq),
        "orientation_mode": "unknown",
    }

    group_reasons = _group_precheck(seq)
    if group_reasons:
        return {"status": "no_solution", **base, "reasons": group_reasons}

    classes, stats = _enumerate_unknown(seq, k, limit=2)

    if not classes:
        if stats.get("limit_truncated"):
            reason = {
                "code": "orientation_search_truncated",
                "message": (
                    "方向未知联合枚举超过搜索状态上限仍未判定，请检查输入批次"
                    "规模或联系维护方"
                ),
                "searched_states": stats["states_visited"],
                "state_limit": SEARCH_STATE_LIMIT,
            }
        else:
            reason = {
                "code": "orientation_constraints_unsatisfiable",
                "message": (
                    "允许每条读数独立取原序列或反向互补、并与闭环次序联合"
                    "穷尽枚举后，不存在每条读数恰好使用一次的闭环；固定方向"
                    "下的度数失衡不能作为本模式的结论（翻转单条读数会同时"
                    "改变两个 k-mer 的入/出度）"
                ),
                "searched_states": stats["states_visited"],
                "pruned": {
                    "no_continuation": stats["no_continuation"],
                    "group_balance": stats["group_balance_pruned"],
                    "endpoint_connectivity": stats[
                        "endpoint_connectivity_pruned"
                    ],
                    "end_gap": stats["end_gap"],
                },
                "dead_end_examples": stats["dead_ends"],
            }
        return {"status": "no_solution", **base, "reasons": [reason]}

    if len(classes) == 1:
        canon, items, rel = classes[0]
        w = _unknown_witness(canon, items, rel, seq, k)
        return {
            "status": "unique",
            **base,
            **w,
            "normalizations": [
                "cyclic_rotation",
                "reverse_complement",
                "per_read_orientation",
            ],
        }

    witnesses = [
        _unknown_witness(canon, items, rel, seq, k)
        for canon, items, rel in classes
    ]
    return {
        "status": "ambiguous",
        **base,
        "message": (
            "即使允许逐条读数独立定向，仍存在多个不同规范等价类的闭环拼法，"
            "给出两条不同规范见证；同一条码的多种等价定向已合并"
        ),
        "witnesses": witnesses,
        "normalizations": [
            "cyclic_rotation",
            "reverse_complement",
            "per_read_orientation",
        ],
    }


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------

def assemble(sequences: list[str], unknown_orientation: bool = False) -> dict:
    """主入口：返回唯一 / 歧义 / 无解三类结果之一。

    unknown_orientation=False（默认）时输入、结论与证据结构与历史版本完全
    兼容；True 时每条读数可独立按原序列或反向互补参与拼接。
    """
    seq = list(sequences)
    k = validate(seq)

    if unknown_orientation:
        return _assemble_unknown(seq, k)

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

    def witness(canon: str, order0: list[int], relation: str) -> dict:
        return {
            "canonical_barcode": canon,
            "barcode": _barcode_from_path(order0, edges, k),
            "canonical_relation": relation,
            "order": [i + 1 for i in order0],
            "evidence": _evidence(order0, edges, k),
        }

    if len(classes) == 1:
        canon, order0, rel = classes[0]
        return {
            "status": "unique",
            "overlap_length": k,
            "read_count": len(edges),
            "canonical_barcode": canon,
            "barcode": _barcode_from_path(order0, edges, k),
            "canonical_relation": rel,
            "order": [i + 1 for i in order0],
            "evidence": _evidence(order0, edges, k),
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

    witnesses = [witness(c, o, r) for c, o, r in classes]
    return {
        "status": "ambiguous",
        "overlap_length": k,
        "read_count": len(edges),
        "message": "存在多个不同规范等价类的闭环拼法，给出两条不同规范见证",
        "witnesses": witnesses,
        "normalizations": ["cyclic_rotation", "reverse_complement"],
    }
