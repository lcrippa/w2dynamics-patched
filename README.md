complex     :  complex QMC code, with ponons and mean-field, working
real        :  real QMC code, with phonons and mean-field, working
experimental:  branches with experimental tweaks
    real_with_dc_rotation: rotating the DC matrix in another basis (for strain)

patches_for_gitlab: these are patches for the real and complex branches of the w2dynamic
master repo. If applied on the given branches, they should include out of the box all mean-field
terms and phonons. 

Correct way of fixing bugs:

- clone the appropriate branch from the w2dynamics gitlab repo
- apply the relative real/complex patch here
- does it apply correctly?
- YES:
    - make your modifications
    - git add -A
    - git diff > NEW_PATCH_NAME.patch
    - save the new patch as real/complex in the appropriate folder
    - git reset
    - clean your repo
    - please port the same changes (by hand) in the complex/real repos here, for safety and consistency
- NO:
    - try to fix conflicts by hand
    - go to YES
