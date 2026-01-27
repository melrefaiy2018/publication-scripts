import numpy as np

def load_hamiltonian(path_file, map_chains_from_saved=None):
    file = open(path_file)
    list_names = file.readline()
    list_names = list_names.strip().split(',')[1:]
    h2_ham = np.loadtxt(path_file, skiprows=1, usecols=np.arange(1, len(list_names)+1), delimiter=',')
    if map_chains_from_saved is not None:
        list_names_new = map_chain_names(list_names, map_chains_from_saved)
    else:
        list_names_new = list_names
    return list_names_new, h2_ham


def map_chain_names(list_names, map_chains_from_saved):
    list_names_new = []
    for name in list_names:
        if map_chains_from_saved is not None:
            if name.split('_')[-3] in map_chains_from_saved.keys():
                list_name_parts = name.split('_')
                list_name_parts[-3] = map_chains_from_saved[list_name_parts[-3]]
                name = '_'.join(list_name_parts)
        list_names_new.append(name)
    return list_names_new


def match_pigment_names(list_name_saved, list_name_site):
    """
    Match pigment names between a saved list and a site list by attempting multiple naming formats.
    
    This function tries to match pigment names by:
    1. Direct exact match
    2. Match with a common prefix extracted from the site list
    
    Args:
        list_name_saved (list): List of pigment names from saved data. 
                                Example: ['CLA_501', 'CLA_502']
        list_name_site (list): List of pigment names from site data.
                               Example: ['PSI_ISIA_A_CLA_501', 'PSI_ISIA_A_CLA_502']
    
    Returns:
        list or None: 
            - If exact match found: returns list_name_saved unchanged
            - If prefix match found: returns list_name_saved with prefix prepended
            - If no match found: returns None
    
    Raises:
        ValueError: If either list is empty
        ValueError: If list_name_site contains strings without at least 3 underscore-separated parts
    
    Example:
        >>> match_pigment_names(['CLA_501'], ['CLA_501', 'CLA_502'])
        ['CLA_501']
        
        >>> match_pigment_names(['CLA_501'], ['PSI_ISIA_A_CLA_501', 'PSI_ISIA_A_CLA_502'])
        ['PSI_ISIA_A_CLA_501']
        
        >>> match_pigment_names(['XYZ_999'], ['PSI_ISIA_A_CLA_501'])
        None
    """
    # Input validation
    if not list_name_saved:
        raise ValueError("list_name_saved cannot be empty")
    if not list_name_site:
        raise ValueError("list_name_site cannot be empty")
    if not isinstance(list_name_saved, (list, tuple)):
        raise TypeError(f"list_name_saved must be a list or tuple, got {type(list_name_saved)}")
    if not isinstance(list_name_site, (list, tuple)):
        raise TypeError(f"list_name_site must be a list or tuple, got {type(list_name_site)}")
    
    # Validate that list_name_site has proper naming format
    first_site_name = list_name_site[0]
    parts = first_site_name.split('_')
    if len(parts) < 3:
        raise ValueError(
            f"Pigment names in list_name_site must have at least 3 underscore-separated parts. "
            f"Got: '{first_site_name}' with {len(parts)} parts"
        )
    
    # Extract common prefix from site names (everything except last 3 parts)
    pre_name = '_'.join(parts[:-3])
    first_saved_name = list_name_saved[0]
    
    # Try to match pigment names
    if first_saved_name in list_name_site:
        # Exact match found
        return list_name_saved
    elif f'{pre_name}_{first_saved_name}' in list_name_site:
        # Prefix match found - prepend prefix to all names in the saved list
        return [f'{pre_name}_{name}' for name in list_name_saved]
    else:
        # No match found
        return None


def update_pigment_name(input_list, replacement_string):
    """
    Processes a list of strings, prepending to the last three parts and replacing others.

    Args:
        input_list: The list of strings to process. Ex: ['PSI_ISIA_A_CLA_501]
        replacement_string: The string to prepend to the last three parts or replace others. Ex: ['protein_A_CLA_501]

    Returns:
        A new list with the processed strings.
    """

    output_list = []

    for item in input_list:
        parts = item.split('_')
        if len(parts) >= 3:
            last_three = parts[-3:]
            new_item = replacement_string + "_" + "_".join(last_three)  # Prepend once to the joined string
            output_list.append(new_item)
        else:
            output_list.append(replacement_string)

    return output_list
