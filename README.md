# supercentenarian
Code accompanying the paper "Human mortality at extreme age"

Three databases are used, `italcent.rda`, `IDL2016.rda` and `francent.rda`

The first can be bought for a small fee from the
National Institute of Statistics by registering at the Contact Center
(https://contact.istat.it) and mentioning the semi-supercentenarian
Survey and Marco Marsili as contact person.


The other two datasets can be obtained freely by registering on http://www.supercentenarians.org/, but the data are continuously updated and may fail to match the version used in the paper. The data version used can be obtained upon request by forwarding the IDL proof of registration to Leo Belzile.

The main script, `Semi_supercentenarian.R`, contains a description of the other files. The file `Semi_supercentenarian_fn.R` contains auxiliary functions (likelihoods, optimisation routine, etc.) called by the user.

