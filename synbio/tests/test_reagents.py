from synbio.reagents import *

class TestReagentRegistry:
    def test_add_by_name(self):
        # define reagents
        r = pure_registry
        r.add_by_name("QD")
        r.add_by_name("test1 test2".split())

class TestMixture:
    def test_mixture_from_recipes(self):
        r = pure_registry
        r.add_by_name("QD")
        # define mixtures
        pos_recipe = dict(PURE.recipe)

        neg_recipe = dict(PURE.recipe)
        neg_recipe.update({r["H20"]: 2, r["DNA"]: 0})

        mixtures = mixtures_from_recipes({
            "pos": pos_recipe,
            "neg": neg_recipe,
            "QD": r["QD"].recipe
        })

