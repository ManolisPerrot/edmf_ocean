#!/bin/bash
# this script run tests to ensure that new developpments 
# do not break previous scripts

# set -e
# set +x

#!/bin/bash

# List of Python scripts
python_scripts=('W005_C500_NO_COR72h_profile_LES_vs_EDMF.py' 'W005_C500_engy_budgets.py' 'WANG1_edmf_vs_les.py' 'WANG1_FR_engy_budgets.py')

#python_scripts=('W005_C500_engy_budgets.py')


# Iterate over the list of Python scripts and execute them
for script in "${python_scripts[@]}"; do
    echo "Testing $script..."
    python "$script"
done


