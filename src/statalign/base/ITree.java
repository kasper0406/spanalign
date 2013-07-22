package statalign.base;

import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.plugins.contree.CNetwork;

public interface ITree {
    double getLogLike();
    double getLogPrior();
    double getOrphanLogLike();
    State getState();

    SubstitutionModel getSubstitutionModel();

    // Methods to get parameters from the TKF92 model.
    double getR();
    double getLambda();
    double getMu();

    double getHeat();
}
