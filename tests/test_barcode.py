"""barcode 核心算法单元测试 + 独立暴力参照交叉验证。"""
import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "app"))

from barcode import (  # noqa: E402
    canonical_circle,
    reverse_complement,
    ValidationError,
    assemble,
)


def reads_of(circle: str, length: int) -> list[str]:
    """从环状串生成每个位置的定长读数（允许 length 大于环长）。"""
    n = len(circle)
    doubled = circle * (length // n + 2)
    return [doubled[i:i + length] for i in range(n)]


def brute_classes(reads: list[str]) -> set[str]:
    """独立参照：枚举边实例的全部排列，收集不同规范环。"""
    n = len(reads)
    classes = set()
    for perm in itertools.permutations(range(n)):
        ok = True
        for a, b in zip(perm, perm[1:]):
            if reads[a][1:] != reads[b][:-1]:
                ok = False
                break
        if ok and reads[perm[-1]][1:] != reads[perm[0]][:-1]:
            ok = False
        if ok:
            circle = reads[perm[0]] + "".join(reads[i][-1] for i in perm[1:])
            classes.add(canonical_circle(circle[:n]))
    return classes


def brute_unknown_classes(reads: list[str]) -> set[str]:
    """方向未知参照：锚定首条正向（整环反向互补等价），枚举其余 2^(n-1)
    种定向与全部 (n-1)! 种次序，收集不同规范环。"""
    n = len(reads)
    rcs = [reverse_complement(s) for s in reads]
    classes = set()
    for bits in range(1 << (n - 1)):
        labels = [reads[0]]
        for j in range(1, n):
            labels.append(reads[j] if (bits >> (j - 1)) & 1 == 0 else rcs[j])
        for perm in itertools.permutations(range(1, n)):
            order = [0] + list(perm)
            if all(
                labels[order[p]][1:] == labels[order[(p + 1) % n]][:-1]
                for p in range(n)
            ):
                circle = labels[0] + "".join(
                    labels[order[i]][-1] for i in range(1, n)
                )
                classes.add(canonical_circle(circle[:n]))
    return classes



def check_witness(test: unittest.TestCase, w: dict, reads: list[str], k: int) -> None:
    """见证内部一致性：次序、重叠证据、条码窗口、规范代表。"""
    n = len(reads)
    length = k + 1
    test.assertEqual(sorted(w["order"]), list(range(1, n + 1)))
    test.assertEqual(len(w["evidence"]), n)
    reps = (n + length) // n + 1  # 读数可能比环长，需要足够多圈
    doubled = w["barcode"] * reps
    for pos, ev in enumerate(w["evidence"]):
        test.assertEqual(ev["position"], pos)
        a = reads[ev["prev"] - 1]
        b = reads[ev["next"] - 1]
        test.assertEqual(ev["prev"], w["order"][pos])
        test.assertEqual(ev["next"], w["order"][(pos + 1) % n])
        test.assertEqual(a[1:], b[:-1])
        test.assertEqual(ev["overlap"], a[1:])
        test.assertEqual(ev["overlap_length"], k)
        test.assertEqual(ev["appended_base"], b[-1])
        test.assertEqual(doubled[pos:pos + length], reads[w["order"][pos] - 1])
    test.assertEqual(canonical_circle(w["barcode"]), w["canonical_barcode"])


def oriented_read(w: dict, reads: list[str], pos: int) -> str:
    """见证第 pos 位实际采用方向后的读数。"""
    s = reads[w["order"][pos] - 1]
    return s if w["orientations"][pos] == "forward" else reverse_complement(s)


def check_unknown_witness(
    test: unittest.TestCase, w: dict, reads: list[str], k: int
) -> None:
    """方向未知见证的逐项一致性：采用方向、定向后序列、相邻重叠。"""
    n = len(reads)
    length = k + 1
    test.assertEqual(sorted(w["order"]), list(range(1, n + 1)))
    test.assertEqual(len(w["evidence"]), n)
    test.assertEqual(len(w["orientations"]), n)
    test.assertEqual(w["canonical_barcode"], w["barcode"])
    reps = (n + length) // n + 1
    doubled = w["barcode"] * reps
    labels = [oriented_read(w, reads, pos) for pos in range(n)]
    for pos, ev in enumerate(w["evidence"]):
        test.assertEqual(ev["position"], pos)
        test.assertEqual(ev["read"], w["order"][pos])
        test.assertEqual(ev["prev"], w["order"][pos])
        test.assertEqual(ev["next"], w["order"][(pos + 1) % n])
        test.assertEqual(ev["orientation"], w["orientations"][pos])
        test.assertEqual(ev["original_sequence"], reads[w["order"][pos] - 1])
        test.assertEqual(ev["oriented_sequence"], labels[pos])
        nxt = labels[(pos + 1) % n]
        test.assertEqual(labels[pos][1:], nxt[:-1])
        test.assertEqual(ev["overlap"], labels[pos][1:])
        test.assertEqual(ev["overlap_length"], k)
        test.assertEqual(ev["appended_base"], nxt[-1])
        test.assertEqual(doubled[pos:pos + length], labels[pos])
    test.assertEqual(canonical_circle(w["barcode"]), w["canonical_barcode"])


class CanonicalTests(unittest.TestCase):
    def test_rotation_and_reverse_complement(self):
        s = "AATCGCAGT"
        for r in range(len(s)):
            rot = s[r:] + s[:r]
            self.assertEqual(canonical_circle(rot), canonical_circle(s))
        self.assertEqual(canonical_circle(reverse_complement(s)), canonical_circle(s))
        self.assertEqual(reverse_complement(reverse_complement(s)), s)
        self.assertEqual(canonical_circle("AACAAC"), "AACAAC")


class ValidationTests(unittest.TestCase):
    def test_count_bounds(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 5)
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 25)

    def test_length_bounds(self):
        with self.assertRaises(ValidationError):
            assemble(["AA"] * 6)
        with self.assertRaises(ValidationError):
            assemble(["AAAAAAAAA"] * 6)

    def test_equal_length_and_alphabet(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA", "AAC", "ACA", "CAA", "ATT", "TT"])  # 不等长
        with self.assertRaises(ValidationError):
            assemble(["AAA", "AAC", "ACA", "CAN", "AAA", "AAT"])  # 非法碱基 N

    def test_lowercase_is_normalized(self):
        r = assemble(["aaa", "aac", "aca", "caa", "aat", "ata"])
        self.assertNotEqual(r["status"], "invalid_input")


class UniqueTests(unittest.TestCase):
    READS = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]  # 环 AATCGC

    def test_unique_basic(self):
        r = assemble(list(self.READS))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        self.assertEqual(r["overlap_length"], 2)
        check_witness(self, r, self.READS, 2)

    def test_reverse_complement_input_same_class(self):
        rc_reads = [reverse_complement(s) for s in self.READS]
        r1 = assemble(list(self.READS))
        r2 = assemble(rc_reads)
        self.assertEqual(r2["status"], "unique")
        self.assertEqual(r1["canonical_barcode"], r2["canonical_barcode"])
        check_witness(self, r2, rc_reads, 2)

    def test_shuffled_input_same_barcode_stable_order(self):
        shuffled = ["TCG", "CAA", "AAT", "GCA", "ATC", "CGC"]
        r = assemble(list(shuffled))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        # 按输入序号稳定分配：起点边 AAT 在打乱后的序号为 3
        self.assertEqual(r["order"][0], 3)
        check_witness(self, r, shuffled, 2)
        # 重复运行结果一致
        self.assertEqual(assemble(list(shuffled))["order"], r["order"])

    def test_duplicate_reads_interchange_not_ambiguous(self):
        # 周期环 AAC：每个读数出现两次，平行边互换只能得到一个规范类
        reads = ["AAC", "ACA", "CAA", "AAC", "ACA", "CAA"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AACAAC")
        self.assertEqual(r["order"], [1, 2, 3, 4, 5, 6])  # 最小可用序号
        check_witness(self, r, reads, 2)

        shuffled = ["ACA", "CAA", "AAC", "AAC", "ACA", "CAA"]
        r2 = assemble(list(shuffled))
        self.assertEqual(r2["status"], "unique")
        self.assertEqual(r2["canonical_barcode"], "AACAAC")
        self.assertEqual(r2["order"], [3, 1, 2, 4, 5, 6])
        check_witness(self, r2, shuffled, 2)

    def test_length_8_small_circle(self):
        # 环长 6 小于读数长度 8：重复结构仍可闭环
        reads = reads_of("AAAAAC", 8)
        self.assertEqual(len(reads), 6)
        r = assemble(list(reads))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AAAAAC")
        check_witness(self, r, reads, 7)


class AmbiguousTests(unittest.TestCase):
    # θ 图：关节点 AA 与 CC 之间各有两条不同路径，去程/回程的交错方式
    # （P+R、Q+S 与 P+S、Q+R）产生不能经旋转或反向互补归一的环。
    #   P: AA→AC→CC   AAC, ACC        R: CC→CA→AA  CCA, CAA
    #   Q: AA→AG→GC→CC AAG, AGC, GCC  S: CC→CG→GA→AA CCG, CGA, GAA
    READS = ["AAC", "ACC", "AAG", "AGC", "GCC",
             "CCA", "CAA", "CCG", "CGA", "GAA"]

    def test_two_distinct_canonical_witnesses(self):
        r = assemble(list(self.READS))
        self.assertEqual(r["status"], "ambiguous")
        self.assertEqual(len(r["witnesses"]), 2)
        barcodes = {w["canonical_barcode"] for w in r["witnesses"]}
        self.assertEqual(len(barcodes), 2)
        for w in r["witnesses"]:
            check_witness(self, w, self.READS, 2)

    def test_matches_brute_force(self):
        expected = brute_classes(self.READS)
        self.assertGreaterEqual(len(expected), 2)
        got = {w["canonical_barcode"] for w in assemble(list(self.READS))["witnesses"]}
        self.assertTrue(got <= expected)


class NoSolutionTests(unittest.TestCase):
    def test_degree_imbalance(self):
        reads = ["AAA", "AAA", "AAC", "ACA", "CAA", "CAA"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertIn("degree_imbalance", codes)
        self.assertNotIn("fragmented_graph", codes)
        imba = next(x for x in r["reasons"] if x["code"] == "degree_imbalance")
        bad = {v["vertex"] for v in imba["vertices"]}
        self.assertEqual(bad, {"AA", "CA"})

    def test_fragmented_graph(self):
        reads = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertIn("fragmented_graph", codes)
        self.assertNotIn("degree_imbalance", codes)
        frag = next(x for x in r["reasons"] if x["code"] == "fragmented_graph")
        self.assertEqual(frag["component_count"], 3)
        self.assertEqual(
            sorted(c[0] for c in frag["components"]), ["AA", "CC", "GG"]
        )
        self.assertTrue(all(len(c) == 1 for c in frag["components"]))

    def test_imbalance_and_fragmentation_reported_together(self):
        reads = ["AAA", "AAA", "AAC", "CCC", "GGG", "TTT"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"degree_imbalance", "fragmented_graph"})
        frag = next(x for x in r["reasons"] if x["code"] == "fragmented_graph")
        # 分量 {AA, AC}（含两个顶点）以及 {CC}、{GG}、{TT}
        self.assertEqual(frag["component_count"], 4)
        self.assertEqual(
            sorted(frag["components"]),
            [["AA", "AC"], ["CC"], ["GG"], ["TT"]],
        )


class BruteForceCrossCheck(unittest.TestCase):
    """随机环状游走生成的小图：求解器结论必须与全排列暴力枚举一致。"""

    def test_random_circles(self):
        rng = random.Random(20260929)
        for _ in range(60):
            n = rng.randint(6, 7)
            circle = "".join(rng.choice("AC") for _ in range(n))
            reads = reads_of(circle, 3)
            expected = brute_classes(reads)
            r = assemble(list(reads))
            if len(expected) == 1:
                self.assertEqual(r["status"], "unique", reads)
                self.assertEqual(r["canonical_barcode"], next(iter(expected)))
                check_witness(self, r, reads, 2)
            else:
                self.assertEqual(r["status"], "ambiguous", reads)
                got = {w["canonical_barcode"] for w in r["witnesses"]}
                self.assertEqual(len(got), 2)
                self.assertTrue(got <= expected)


class UnknownOrientationTests(unittest.TestCase):
    """方向未知模式：每条读数可独立正向或反向互补，定向与次序联合求解。"""

    READS = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]  # 环 AATCGC

    def test_forward_only_submission(self):
        r = assemble(list(self.READS), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        self.assertEqual(r["orientation_mode"], "unknown")
        self.assertEqual(r["orientations"], ["forward"] * 6)
        check_unknown_witness(self, r, self.READS, 2)

    def test_per_read_flips_recovered_jointly(self):
        # 随机翻转其中 4 条：不能靠整批翻转，必须逐条定向
        flipped = [reverse_complement(s) if i in (0, 3, 4, 5) else s
                   for i, s in enumerate(self.READS)]
        r = assemble(list(flipped), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        self.assertEqual(sorted(r["order"]), [1, 2, 3, 4, 5, 6])
        check_unknown_witness(self, r, flipped, 2)
        # 被翻转的读数在见证中必须以 reverse_complement 恢复回环方向
        by_read = {r["order"][p]: r["orientations"][p]
                   for p, _ in enumerate(r["order"])}
        for i in (0, 3, 4, 5):
            self.assertEqual(by_read[i + 1], "reverse_complement")
        self.assertEqual(by_read[2], "forward")

    def test_single_flip_that_fixed_mode_cannot_close(self):
        # 只翻转一条读数：固定方向报度数失衡，方向未知可闭环
        broken = list(self.READS)
        broken[0] = reverse_complement(broken[0])  # AAT -> ATT
        rf = assemble(list(broken))
        self.assertEqual(rf["status"], "no_solution")
        self.assertIn("degree_imbalance",
                      {x["code"] for x in rf["reasons"]})
        ru = assemble(list(broken), unknown_orientation=True)
        self.assertEqual(ru["status"], "unique")
        self.assertEqual(ru["canonical_barcode"], "AATCGC")
        check_unknown_witness(self, ru, broken, 2)

    def test_each_index_used_exactly_once(self):
        flipped = [reverse_complement(s) for s in self.READS]
        r = assemble(list(flipped), unknown_orientation=True)
        self.assertEqual(sorted(r["order"]), list(range(1, 7)))
        # 即便全部反向互补，也是逐条定向的结果，而非整批翻转的特判
        self.assertEqual(set(r["orientations"]), {"reverse_complement"})
        check_unknown_witness(self, r, flipped, 2)

    def test_palindromic_reads_direction_collapsed(self):
        # AAA 回文：正反向序列相同，见证统一记 forward，且可拼 AAAAAA
        reads = ["AAA"] * 6
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AAAAAA")
        self.assertEqual(r["orientations"], ["forward"] * 6)
        check_unknown_witness(self, r, reads, 2)

    def test_equivalent_orientations_not_ambiguous(self):
        # 周期环 ATATAT：每条读数的两个方向都能参与同一闭环，
        # 同一条码的多种等价定向不得制造歧义。
        reads = reads_of("ATATAT", 3)
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "ATATAT")
        check_unknown_witness(self, r, reads, 2)

    def test_duplicate_reads_with_mixed_flips(self):
        # 周期环 AAC 的两组重复读数，分别做逐条翻转
        reads = ["AAC", "ACA", "CAA", "AAC", "ACA", "CAA"]
        flipped = [reverse_complement(s) if i in (0, 2, 4) else s
                   for i, s in enumerate(reads)]
        r = assemble(list(flipped), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AACAAC")
        check_unknown_witness(self, r, flipped, 2)

    def test_ambiguous_still_two_stable_witnesses(self):
        # 即便允许逐条定向，两个规范类仍然不同
        reads = ['ATT', 'AGA', 'TTA', 'TAA', 'GAT', 'TAT', 'TAT', 'AAG']
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "ambiguous")
        self.assertEqual(len(r["witnesses"]), 2)
        barcodes = {w["canonical_barcode"] for w in r["witnesses"]}
        self.assertEqual(barcodes, {"AAGATAAT", "AAGATATT"})
        for w in r["witnesses"]:
            check_unknown_witness(self, w, reads, 2)
        # 稳定：重复调用结果一致
        r2 = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(
            [w["order"] for w in r2["witnesses"]],
            [w["order"] for w in r["witnesses"]],
        )

    def test_group_odd_degree_reason(self):
        # AAT/TAA 的配对组 {AA,TT} 度 3（奇），{AT} 度 3（奇）
        reads = ["AAT", "TAA", "AAT", "TAA", "AAT", "TAA"]
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"orientation_group_imbalance"})
        odd = r["reasons"][0]["odd_degree_groups"]
        groups = {tuple(g["group"]): g["degree"] for g in odd}
        self.assertEqual(groups, {("AT",): 3, ("TA",): 3})

    def test_group_fragmentation_reason(self):
        reads = ["AAA", "AAA", "AAA", "CCC", "CCC", "CCC"]
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"orientation_group_fragmented"})
        self.assertEqual(r["reasons"][0]["component_count"], 2)

    def test_exhaustive_unsatisfiable_not_degree_imbalance(self):
        # 通过组级必要条件、但任何逐条定向都无法闭环：
        # 必须给方向约束不可满足的可复核原因，而不是固定方向度数失衡。
        reads = ['GTC', 'TCT', 'ACG', 'GCG', 'GCC', 'AGG']
        rf = assemble(list(reads))
        self.assertEqual(rf["status"], "no_solution")  # 固定方向确有度数失衡
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"orientation_constraints_unsatisfiable"})
        self.assertNotIn("degree_imbalance", codes)
        reason = r["reasons"][0]
        self.assertGreaterEqual(reason["searched_states"], 1)
        self.assertTrue(reason["dead_end_examples"])
        de = reason["dead_end_examples"][0]
        self.assertIn("at_kmer", de)
        self.assertIn("candidate_edges", de)

    def test_shuffled_flips_stable_and_deterministic(self):
        rng = random.Random(20260930)
        for _ in range(20):
            n = rng.randint(6, 9)
            circle = "".join(rng.choice("ACGT") for _ in range(n))
            reads = reads_of(circle, 3)
            reads = [reverse_complement(s) if rng.random() < 0.5 else s
                     for s in reads]
            rng.shuffle(reads)
            r1 = assemble(list(reads), unknown_orientation=True)
            r2 = assemble(list(reads), unknown_orientation=True)
            self.assertEqual(r1["status"], r2["status"])
            if r1["status"] == "unique":
                self.assertEqual(r1["order"], r2["order"])
                self.assertEqual(r1["orientations"], r2["orientations"])
                check_unknown_witness(self, r1, reads, 2)
            elif r1["status"] == "ambiguous":
                for w1, w2 in zip(r1["witnesses"], r2["witnesses"]):
                    self.assertEqual(w1["order"], w2["order"])
                    check_unknown_witness(self, w1, reads, 2)

    def test_brute_force_cross_check(self):
        """随机环+逐条翻转：求解器结论与 2^(n-1)*(n-1)! 暴力枚举一致。"""
        rng = random.Random(20260931)
        for _ in range(40):
            n = rng.randint(6, 7)
            circle = "".join(rng.choice("AC") for _ in range(n))
            reads = reads_of(circle, 3)
            reads = [reverse_complement(s) if rng.random() < 0.5 else s
                     for s in reads]
            rng.shuffle(reads)
            expected = brute_unknown_classes(reads)
            self.assertTrue(expected)
            r = assemble(list(reads), unknown_orientation=True)
            if len(expected) == 1:
                self.assertEqual(r["status"], "unique", reads)
                self.assertEqual(r["canonical_barcode"], next(iter(expected)))
                check_unknown_witness(self, r, reads, 2)
            else:
                self.assertEqual(r["status"], "ambiguous", reads)
                got = {w["canonical_barcode"] for w in r["witnesses"]}
                self.assertEqual(len(got), 2)
                self.assertTrue(got <= expected)


class CompatibilityTests(unittest.TestCase):
    """未启用方向未知时，原有输入、结论与证据结构完全兼容。"""

    READS = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]

    def test_default_equals_false(self):
        a = assemble(list(self.READS))
        b = assemble(list(self.READS), unknown_orientation=False)
        self.assertEqual(a, b)
        self.assertNotIn("orientation_mode", a)
        self.assertNotIn("orientations", a)
        # 固定模式证据不含方向字段
        self.assertNotIn("orientation", a["evidence"][0])
        self.assertNotIn("oriented_sequence", a["evidence"][0])

    def test_existing_fixtures_unchanged(self):
        r = assemble(list(self.READS))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        self.assertEqual(r["order"], [1, 2, 3, 4, 5, 6])
        check_witness(self, r, self.READS, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
