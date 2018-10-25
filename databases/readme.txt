-Genes.csv:
genes name dictionary to match up the entrez gene id

-SBSM_seperated.csv:
CHEMBLID, DrugBankID, LigExpoID, SBSMID...etc for each compound

-comp_names.pickle:
A dictionary to match the SBSMID to the drug's name or CHEMBLID or DrugBankID depending on which is available.

-final_database.csv:
this file contains all the compounds in our database in the following format;

SBSMID	SMILES	TARGETS

-final_database.pickle:
created from final_database.csv for faster loading into hitpick.py
when loaded, it is unpickled into a dictionary:
{  SBSMID : [ SMILES, INCHI, bitvector, targets]   }

-precision_table.pickle:
a pickled dictionary storing the precision score of a targets based on rank, Tanimoto sim and # of occur.

-bayesian_models:
all the models of the 2738 genes

in case of any question:
email: sabrihamad@outlook.com
