import numpy as np
from colorama import Fore, Style

def is_in_range(value, value_name, reference, tolerance):
    if np.abs(value-reference) > tolerance:
        print(Fore.RED+'WARNING: test failed!'+value_name,'is', value, ', but the reference value is', reference,'+-', tolerance)
        print(Style.RESET_ALL)
    else:
        print(Fore.GREEN+'Test sucessfully passed',value_name,'is', value, ', and the reference value is', reference,'+-', tolerance)
        print(Style.RESET_ALL)