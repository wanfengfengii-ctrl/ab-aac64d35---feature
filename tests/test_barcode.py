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


def brute_classes_unknown(reads: list[str]) -> set[str]:
    """方向未知参照：枚举 2^n 逐条定向 × 边实例全排列，收集不同规范环。"""
    n = len(reads)
    classes = set()
    for bits in range(1 << n):
        rs = [
            reverse_complement(reads[i]) if (bits >> i) & 1 else reads[i]
            for i in range(n)
        ]
        for perm in itertools.permutations(range(n)):
            if all(
                rs[perm[i]][1:] == rs[perm[(i + 1) % n]][:-1]
                for i in range(n)
            ):
                circle = rs[perm[0]] + "".join(
                    rs[perm[j]][-1] for j in range(1, n)
                )
                classes.add(canonical_circle(circle[:n]))
    return classes


def check_unknown_witness(
    test: unittest.TestCase, w: dict, reads: list[str], k: int
) -> None:
    """方向未知见证一致性：方向、定向后序列、闭环重叠、每个序号一次。"""
    n = len(reads)
    test.assertEqual(sorted(w["order"]), list(range(1, n + 1)))
    test.assertEqual(len(w["evidence"]), n)
    test.assertEqual(len(w["read_orientations"]), n)
    oriented = {}
    for ro in w["read_orientations"]:
        i = ro["read"] - 1
        test.assertIn(ro["orientation"], ("forward", "reverse_complement"))
        expected = (
            reads[i]
            if ro["orientation"] == "forward"
            else reverse_complement(reads[i])
        )
        test.assertEqual(ro["oriented_sequence"], expected)
        oriented[ro["read"]] = expected
    for pos, ev in enumerate(w["evidence"]):
        test.assertEqual(ev["position"], pos)
        a, b = oriented[ev["prev"]], oriented[ev["next"]]
        test.assertEqual(ev["oriented_sequence"], a)
        test.assertEqual(a[1:], b[:-1])
        test.assertEqual(ev["overlap"], a[1:])
        test.assertEqual(ev["overlap_length"], k)
    test.assertEqual(canonical_circle(w["barcode"]), w["canonical_barcode"])


class UnknownOrientationTests(unittest.TestCase):
    READS = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]  # 环 AATCGC

    def test_unknown_flag_off_is_unchanged(self):
        r = assemble(list(self.READS), unknown_orientation=False)
        self.assertEqual(r["status"], "unique")
        self.assertNotIn("orientation_mode", r)
        self.assertNotIn("orientation", r["evidence"][0])
        self.assertNotIn("read_orientations", r)

    def test_already_forward_batch(self):
        r = assemble(list(self.READS), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["orientation_mode"], "unknown")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        self.assertEqual(
            [ro["orientation"] for ro in r["read_orientations"]],
            ["forward"] * 6,
        )
        check_unknown_witness(self, r, self.READS, 2)

    def test_randomly_rc_flipped_reads_recover_same_barcode(self):
        rng = random.Random(20260930)
        for _ in range(20):
            flipped = [
                reverse_complement(s) if rng.random() < 0.5 else s
                for s in self.READS
            ]
            r = assemble(list(flipped), unknown_orientation=True)
            self.assertEqual(r["status"], "unique", flipped)
            self.assertEqual(r["canonical_barcode"], "AATCGC", flipped)
            check_unknown_witness(self, r, flipped, 2)

    def test_global_flip_is_one_class_not_ambiguous(self):
        # 整批反向互补必须仍只产生一个规范类
        flipped = [reverse_complement(s) for s in self.READS]
        r = assemble(list(flipped), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")

    def test_two_classes_return_two_stable_witnesses(self):
        theta = ["AAC", "ACC", "AAG", "AGC", "GCC",
                 "CCA", "CAA", "CCG", "CGA", "GAA"]
        rng = random.Random(7)
        flipped = [
            reverse_complement(s) if rng.random() < 0.5 else s
            for s in theta
        ]
        r = assemble(list(flipped), unknown_orientation=True)
        self.assertEqual(r["status"], "ambiguous")
        self.assertEqual(len(r["witnesses"]), 2)
        barcodes = {w["canonical_barcode"] for w in r["witnesses"]}
        self.assertEqual(len(barcodes), 2)
        for w in r["witnesses"]:
            check_unknown_witness(self, w, flipped, 2)
        # 稳定：重复运行两份见证完全一致
        r2 = assemble(list(flipped), unknown_orientation=True)
        self.assertEqual(
            [w["canonical_barcode"] for w in r["witnesses"]],
            [w["canonical_barcode"] for w in r2["witnesses"]],
        )
        self.assertEqual(
            [w["order"] for w in r["witnesses"]],
            [w["order"] for w in r2["witnesses"]],
        )

    def test_unsatisfiable_gives_orientation_reason_not_degree(self):
        # 固定方向下是 degree_imbalance；方向未知下必须给出方向约束不可满足
        reads = ["AAA", "AAA", "AAC", "ACA", "CAA", "CAA"]
        known = assemble(list(reads))
        self.assertEqual(known["status"], "no_solution")
        self.assertIn(
            "degree_imbalance", {x["code"] for x in known["reasons"]}
        )
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"orientation_constraints_unsatisfiable"})
        self.assertNotIn("degree_imbalance", codes)
        reason = r["reasons"][0]
        self.assertIn("necessary_conditions", reason)

    def test_folded_fragmentation_is_orientation_independent(self):
        reads = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "no_solution")
        self.assertEqual(
            r["reasons"][0]["code"], "orientation_fragmented_graph"
        )

    def test_palindromic_reads(self):
        # 偶数长度回文读数：定向不改变图
        reads = ["TATA", "ATAT", "TATA", "ATAT", "TATA", "ATAT"]
        r = assemble(list(reads), unknown_orientation=True)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "ATATAT")
        check_unknown_witness(self, r, reads, 3)


class UnknownBruteForceCrossCheck(unittest.TestCase):
    """2^n 定向 × 全排列暴力参照交叉验证（小图）。"""

    def test_random(self):
        rng = random.Random(20260930)
        for _ in range(40):
            n = 6
            length = rng.choice([3, 4])
            if rng.random() < 0.6:
                circle = "".join(rng.choice("ACGT") for _ in range(n))
                reads = reads_of(circle, length)
            else:
                reads = [
                    "".join(rng.choice("ACGT") for _ in range(length))
                    for _ in range(n)
                ]
            reads = [
                reverse_complement(s) if rng.random() < 0.5 else s
                for s in reads
            ]
            expected = brute_classes_unknown(reads)
            r = assemble(list(reads), unknown_orientation=True)
            if not expected:
                self.assertEqual(r["status"], "no_solution", reads)
            elif len(expected) == 1:
                self.assertEqual(r["status"], "unique", reads)
                self.assertEqual(
                    r["canonical_barcode"], next(iter(expected)), reads
                )
                check_unknown_witness(self, r, reads, length - 1)
            else:
                self.assertEqual(r["status"], "ambiguous", reads)
                got = {w["canonical_barcode"] for w in r["witnesses"]}
                self.assertEqual(len(got), 2)
                self.assertTrue(got <= expected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
