from synbio.experiment import *
from synbio.reagents import *


class TestExperiment:
    def test_add_conditions(self):
        # define reagents
        r = PURE_reagents

        # define mixtures
        pos_recipe = dict(PURE.recipe)

        neg_recipe = dict(PURE.recipe)
        neg_recipe.update({r["H20"]: 2, r["DNA"]: 0})

        mixtures = mixture_from_recipes({
            "pos": pos_recipe,
            "neg": neg_recipe,
        })

        # make experiment object
        exp = Experiment("test_experiment", r)

        # define experimental conditions
        exp.add_conditions(mixtures, replicates=3)  # only defining one
        # optional parameter
