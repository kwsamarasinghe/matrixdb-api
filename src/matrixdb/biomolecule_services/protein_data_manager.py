from itertools import groupby


class ProteinDataManager:

    def __init__(self, database_connection, meta_data_cache):
        self.meta_data_cache = meta_data_cache
        self.database_connection = database_connection

    def get_gene_expressions(self, protein_id):
        expressions_by_protein = self.database_connection["geneExpression"].find_one({
            "uniprot": protein_id
        })

        gene_expressions = list()
        if expressions_by_protein is not None:

            expression_values = list()
            for expression in expressions_by_protein["expressions"]:
                for tissue_id, expression_value in expression.items():
                    expression_by_tissue = expression_value[0]
                    expression_by_tissue['tissueId'] = tissue_id
                    expression_values.append(expression_by_tissue)

            # filtered_expression = filter(lambda e: "life cycle" in e["sex"], expressions_by_protein["expressions"])
            grouped_expressions = groupby(sorted(list(expression_values), key=lambda x: x["uberonName"]),
                                          key=lambda x: x["uberonName"])

            for tissue_uberon_name, group in grouped_expressions:
                max_group = max(group, key=lambda x: x["tpm"])
                gene_expressions.append({
                    "tissueId": max_group["tissueId"],
                    "tissue": tissue_uberon_name.strip('"'),
                    "tpm": max_group["tpm"]
                })
        return gene_expressions

    def get_gene_expression_for_proteins(self, protein_ids):
        expressions_by_proteins = self.database_connection["geneExpression"].find({
            "uniprot": {
                '$in': protein_ids
            }
        })

        gene_expressions = dict()
        if expressions_by_proteins is not None:
            for expressions_by_protein in expressions_by_proteins:
                protein = expressions_by_protein['uniprot']
                gene_expressions[protein] = list()
                expression_values = list()
                for expression in expressions_by_protein["expressions"]:
                    for tissue_id, expression_value in expression.items():
                        expression_by_tissue = expression_value[0]
                        expression_by_tissue['tissueId'] = tissue_id
                        expression_values.append(expression_by_tissue)

                # filtered_expression = filter(lambda e: "life cycle" in e["sex"], expressions_by_protein["expressions"])
                grouped_expressions = groupby(sorted(list(expression_values), key=lambda x: x["uberonName"]),
                                              key=lambda x: x["uberonName"])

                for tissue_uberon_name, group in grouped_expressions:
                    max_group = max(group, key=lambda x: x["tpm"])
                    gene_expressions[protein].append({
                        "tissueId": max_group["tissueId"],
                        "tissue": tissue_uberon_name.strip('"'),
                        "tpm": max_group["tpm"]
                    })
        return gene_expressions

    def get_proteomics_expressions(self, protein_id):

        proteomics_expression_by_protein = self.database_connection["proteomicsExpression"].find_one({
            "uniprot": protein_id
        })

        prot_expressions = list()
        if proteomics_expression_by_protein is not None:
            for e in proteomics_expression_by_protein["expressions"]:

                if "BTO" in e["tissueId"]:
                    if e['tissueId'] in self.meta_data_cache["bto"]:
                        tissueName = self.meta_data_cache["bto"][e["tissueId"]]["name"]

                if "UBERON" in e["tissueId"]:
                    if e['tissueId'] in self.meta_data_cache["uberon"]:
                        tissueName = self.meta_data_cache["uberon"][e["tissueId"]]["name"]

                if tissueName is not None:
                    prot_expressions.append(
                        {
                            "tissue": tissueName,
                            "sampleName": e["sampleName"] if "sampleName" in e else "",
                            "sample": e["sample"] if "sample" in e else "",
                            "score": e["confidenceScore"] if "confidenceScore" in e else 0,
                        })

        return prot_expressions

    def get_proteomics_expressions_for_proteins(self, protein_ids):

        proteomics_expression_by_proteins = self.database_connection["proteomicsExpression"].find({
            "uniprot": {
                '$in': protein_ids
            }
        })

        prot_expressions = dict()
        if proteomics_expression_by_proteins is not None:
            for proteomics_expression_by_protein in list(proteomics_expression_by_proteins):
                prot_expressions[proteomics_expression_by_protein['uniprot']] = list()
                for e in proteomics_expression_by_protein["expressions"]:

                    if "BTO" in e["tissueId"]:
                        if e['tissueId'] in self.meta_data_cache["bto"]:
                            tissueName = self.meta_data_cache["bto"][e["tissueId"]]["name"]

                    if "UBERON" in e["tissueId"]:
                        if e['tissueId'] in self.meta_data_cache["uberon"]:
                            tissueName = self.meta_data_cache["uberon"][e["tissueId"]]["name"]

                    if tissueName is not None:
                        prot_expressions[proteomics_expression_by_protein['uniprot']].append(
                            {
                                "tissue": tissueName,
                                "sampleName": e["sampleName"] if "sampleName" in e else "",
                                "sample": e["sample"] if "sample" in e else "",
                                "score": e["confidenceScore"] if "confidenceScore" in e else 0,
                            })

        return prot_expressions