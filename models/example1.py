# Initialize the constructor
from MANAscore import MANAscore

mana = MANAscore()
mana.load_voting_models()

## predict MANAscore for your own data
mana.predict_and_save_both(
    input_file_i='./3gene/validation/MD01-004_known.csv',
    input_file_ni='./3gene/validation/MD01-004_known_RNA.csv',
    output_file_i='MD01-004_voting_i_predict_score.csv',
    output_file_ni='MD01-004_voting_ni_predict_score.csv'
)

