# Initialize the constructor
from MANAscore import MANAscore

mana = MANAscore()
#mana.load_voting_models()
mana.load_and_split_data()
#mana.train_logistic_models()
mana.train_logistic_models()
#mana.train_random_forest_models()
mana.train_random_forest_models()
## voting models
mana.create_and_fit_voting_classifiers()
## save models
mana.save_voting_models()

## evaluate models
mana.evaluate_classifiers()

## predict MANAscore
mana.predict_and_save_both(
    input_file_i='./3gene/validation/MD01-004_known.csv',
    input_file_ni='./3gene/validation/MD01-004_known_RNA.csv',
    output_file_i='MD01-004_voting_i_predict_score.csv',
    output_file_ni='MD01-004_voting_ni_predict_score.csv'
)

