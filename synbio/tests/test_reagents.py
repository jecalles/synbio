from synbio.reagents import *

class TestMixture:
    def test_mixture_from_recipes(self):
        r = PURE_reagents
        # define mixtures
        pos_recipe = dict(PURE.recipe)

        neg_recipe = dict(PURE.recipe)
        neg_recipe.update({r["H20"]: 2, r["DNA"]: 0})

        mixtures = mixture_from_recipes({
            "pos": pos_recipe,
            "neg": neg_recipe,
        })

