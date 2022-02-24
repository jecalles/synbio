from synbio.reagents import *
from synbio.experiment import *

class TestExperiment:
    def test_add_conditions(self):
        # define reagents
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

        # make experiment object
        exp = Experiment("test_experiment", r)

        # define experimental conditions
        exp.add_conditions(mixtures, replicates=3) # only defining one
        # optional parameter
