import urllib.parse

def convert_attributes(attr_str):
    """
    Convert a GFF attribute string into a Python dictionary.
    
    Parameters:
        attr_str (str): A semicolon-separated attribute string.
        
    Returns:
        dict: A dictionary with keys and their corresponding values.
              If a value contains commas, it is split into a list.
    """
    attr_dict = {}
    # Split the string by semicolon
    for pair in attr_str.split(';'):
        pair = pair.strip()
        if not pair:
            continue
        if '=' in pair:
            key, value = pair.split('=', 1)
            key = key.strip()
            # Decode URL-encoded parts of the value
            value = urllib.parse.unquote(value.strip())
            # If the value contains commas, split it into a list
            if ',' in value:
                value = [v.strip() for v in value.split(',')]
            attr_dict[key] = value
        else:
            # If no '=' is present, store the key with a True flag
            attr_dict[pair] = True
    return attr_dict

# Example usage:
if __name__ == '__main__':
    attribute_string = ("ID=FBgn0031802_d7916e9908;Name=ppk7;"
                        "to_species=Drosophila melanogaster;to_name=ppk12;"
                        "Dbxref=FlyBase:FBgn0034730,FlyBase:FBan0010972,FlyBase_Annotation_IDs:CG10972,"
                        "GB_protein:AAF46849,GB:AY226543,GB_protein:AAO47369,UniProt/TrEMBL:Q9W250,"
                        "INTERPRO:IPR001873,EntrezGene:37566,FlyMine:FBgn0034730,OrthoDB9_1_Diptera:EOG091502CU,"
                        "OrthoDB9_1_Insecta:EOG090W025S,OrthoDB9_1_Arthropoda:EOG090X0232,"
                        "OrthoDB9_1_Metazoa:EOG091G03EH,OrthoDB9_1_Drosophila:EOG091908MQ,"
                        "UniProt/GCRP:Q9W250,KEGG_GENES:dme:Dmel_CG10972,AlphaFold_DB:Q9W250,"
                        "DRscDB:37566/tissue%3DAll,EMBL-EBI_Single_Cell_Expression_Atlas:FBgn0034730,"
                        "MARRVEL_MODEL:37566,FlyAtlas2:FBgn0034730;"
                        "Target=2R 22427479 22429584 +;"
                        "diopt_source=eggNOG,Panther,OrthoInspector,OrthoDB,Phylome,Domainoid")
    
    converted = convert_attributes(attribute_string)
    print(converted['Target'])