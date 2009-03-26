function [test_perp,unseen_perp,lam] = plsa_delinterp(pwgz,pdgz,pz,train,heldout,test,unseen_test)

[ho_lik] = plsa_getlik(pwgz,pdgz,pz,heldout);
[lam,hist,perp] = delinterp(train,heldout,ho_lik',0);

[tst_lik] = plsa_getlik(pwgz,pdgz,pz,test);
[a,b,test_perp] = delinterp(train,test,tst_lik',lam);

[unseen_lik] = plsa_getlik(pwgz,pdgz,pz,unseen_test);
[a,b,unseen_perp] = delinterp(train,unseen_test,unseen_lik',lam);

test_perp = full(test_perp);
unseen_perp = full(unseen_perp);

lam = full(lam);
