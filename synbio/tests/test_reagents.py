from synbio.reagents import *
import synbio.tests.utils as testutils

class TestReagent:
    pure_reagent_names = "A B C D".split()
    pure_reagents = {
        name: Reagent(name)
        for name in pure_reagent_names
    }
    mix1 = Mixture("mix1", {
        pure_reagents["A"]:1,
        pure_reagents["B"]:2,
    })
    mix2 = Mixture("mix2", {
        pure_reagents["C"]:2,
        pure_reagents["D"]:3
    })
    mix3 = Mixture("mix3", {
        mix1:1,
        mix2:2,
        pure_reagents["A"]:3
    })

    registry = mix1._registry

    def test___eq__(self):
        reg = self.registry
        # test simple equality
        a = self.pure_reagents["A"]
        assert a == Reagent("A")

        mix1 = self.mix1
        assert mix1 == reg["mix1"]


class TestRecipe:
    def test___add__(self):
        pure = reagent_registry["PURE"]
        dna = reagent_registry["DNA"]

        QD_1in10 = Recipe({
            Reagent("QD"): 1,   # test Reagent init in Recipe definition
            "H2O": 9            # test str key access for existing reagent
        })
        added = pure.recipe + QD_1in10
        assert added == Recipe({
            'Sol-A': 4,
            'Sol-B': 3,
            dna: 1,             # check reagent input
            'RNase-Inh': 1,
            'H2O': 10,
            'QD': 1
        })
        assert added == pure.recipe + {"QD":1, "H2O": 9} # check dict add

    def test___sub__(self):
        pure: Mixture = reagent_registry["PURE"]
        PURE_dtRNA = pure.recipe - {"Sol-A": 2} + {Reagent("AA"):1, Reagent(
            "tRNA"): 1}
        assert PURE_dtRNA == Recipe({
            'Sol-A': 2,
            'Sol-B': 3,
            'AA':1,
            'tRNA':1,
            'DNA': 1,
            'RNase-Inh': 1,
            'H2O': 1,
        })
        assert pure.recipe - {"H2O": 2} == Recipe({
            'Sol-A': 4,
            'Sol-B': 3,
            'DNA': 1,
            'RNase-Inh': 1,
        })

    def test___mul__(self):
        pure = reagent_registry["PURE"]
        assert pure.recipe * 2.5 == Recipe({
            'Sol-A': 10,
            'Sol-B': 7.5,
            'DNA': 2.5,
            'H2O': 2.5,
            'RNase-Inh': 2.5
        })

    def test___div__(self):
        pure = reagent_registry["PURE"]
        assert pure.recipe / 100 == Recipe({
            'Sol-A': 0.04,
            'Sol-B': 0.03,
            'DNA': 0.01,
            'H2O': 0.01,
            'RNase-Inh': 0.01
        })

    def test_normalized(self):
        pure = reagent_registry["PURE"]
        assert pure.recipe.normalized == Recipe({
            'Sol-A': 0.4,
            'Sol-B': 0.3,
            'DNA': 0.1,
            'H2O': 0.1,
            'RNase-Inh': 0.1
        })

    def test_flattened(self):
        new_rgts = {
            name: Reagent(name)
            for name in "A B C D".split()
        }
        test_recipe = Recipe({
            "A": 2,
            "B": 3,
            Mixture("C+D", {"C": 2, "D": 1}): 3, # should not need normalization
            Mixture("A+'C+D'", {"A": 1, "C+D":1}): 2,
            "D":1
        })
        assert test_recipe.flat == Recipe({
            "A": 2 + 1,
            "B": 3,
            "C": (2/3)*3 + (2/3)*(1/2)*2,
            "D": (1/3)*3 + (1/3)*(1/2)*2 + 1
        })

    def test_get_reagents(self):
        new_rgts = {
            name: Reagent(name)
            for name in "A B C D".split()
        }
        test_mix = Mixture("test_mix",{
            "A": 2,
            "B": 3,
            Mixture("C+D", {"C": 2, "D": 1}): 3,
            # should not need normalization
            Mixture("A+'C+D'", {"A": 1, "C+D": 1}): 2,
            "D": 1
        })
        pure = reagent_registry["PURE"]
        assert flatten_reagents([test_mix, pure]) == {
            reagent_registry["Sol-A"],
            reagent_registry["Sol-B"],
            reagent_registry["H2O"],
            reagent_registry["DNA"],
            reagent_registry["RNase-Inh"],
            reagent_registry["A"],
            reagent_registry["B"],
            reagent_registry["C"],
            reagent_registry["D"],
        }

class TestMixture:
    pass

