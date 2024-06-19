In this github is deposited the code of the windowing method for the frustration analysis of molecular dynamics.

## A contact-based analysis of local energetic frustration dynamics identifies key residues enabling RfaH fold-switch 

![Workflow of the method]([url completa de la imagen]([https://drive.google.com/file/d/1dZJsgxQgcfhFQem1M0hsxkDIxF_E2xzl/view?usp=sharing](https://github.com/pb3lab/RfaH-frustration/blob/main/Figures/diagram.png))

### ABSTRACT
Fold-switching enables metamorphic proteins to reversibly interconvert between two highly dissimilar native states to regulate their protein functions. While about 100 proteins have been identified to undergo fold-switching, unveiling the key residues behind this mechanism for each protein remains challenging. Reasoning that fold-switching in proteins is driven by dynamic changes in local energetic frustration, we combined fold-switching simulations generated using simplified structure-based models (SBMs) with frustration analysis, to identify key residues involved in this process based on the change in the density of minimally frustrated contacts during refolding. Using this approach to analyze the fold-switch of the bacterial transcription factor RfaH, we identified 20 residues that significantly change their frustration during its fold-switch, some of which have been experimentally and computationally reported in previous works. Our approach, which we developed as an additional tool for the frustratometeR package, highlights the role of local frustration dynamics in protein fold-switching and offers a robust tool to enhance our understanding of other proteins with significant conformational shifts.

## Dependencies
- FrustraR R Package
- Python librabries:
    - numpy
    - pandas

### How to run the pipeline
python3 pipeline_for_frustration_clean.py patho_to_pdbs path_to_results #frames #protein_length #sim_start
Example: python3 pipeline_for_frustration_clean.py /home/titanx1/Documents/sims/frustra/rfah/traj1/ /home/titanx1/Documents/sims/frustra/rfah/traj1/Results/ 162 0

#sys.argv[1] -> path to pdbs files
#sys.argv[2] -> path where do you want to save the frustration results
#sys.argv[3] -> numbers o frames
#sys.argv[4] -> start frame
