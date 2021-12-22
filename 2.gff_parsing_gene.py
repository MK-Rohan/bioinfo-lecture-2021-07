#! /data/home/Minkyu/anaconda3/bin/python

import sys

inp = sys.argv[1]
oup = sys.argv[2]

data_lst = list()
product_lst = list()
trans_lst = list()
protein_lst = list()
ID_dic = dict()


def Usage():
    print(
        "python 2.gff_parsing_gene.py <<input file path + name>> <<output file Path>>"
    )


def Parsing(file):
    global data_lst, ID_dic, product_lst, trans_lst, protein_lst
    infile = open(file, "r")
    outfile = open(oup + "gene_level_" + file[file.rfind("/") + 1 :], "w")
    header = "ID\tGeneID\tCGNC\tgbkey\tgene_biotype\tgene\tgenes_synonym\tproduct\ttranscript_id\tprotein_id"
    # print(header)
    outfile.write(header)
    for line in infile:
        if not line.startswith("#"):
            Ln = line.strip().split("\t")
            Ln2 = Ln[8].strip().split(";")[0]
            if Ln2.startswith("ID=gene"):
                if not data_lst:  # infor is empty list
                    GENE(Ln[8])
                else:
                    if len(product_lst) > 1:
                        Product = ("; ").join(product_lst)
                    elif not product_lst:
                        Product = "."
                    else:
                        Product = product_lst[0]
                    if len(trans_lst) > 1:
                        Trans = ("; ").join(trans_lst)
                    elif not trans_lst:
                        Trans = "."
                    else:
                        Trans = trans_lst[0]
                    if len(protein_lst) > 1:
                        Protein = ("; ").join(protein_lst)
                    elif not protein_lst:
                        Protein = "."
                    else:
                        Protein = protein_lst[0]
                    data_lst += Product, Trans, Protein
                    out_data = (("&\t&").join(data_lst)).strip().split("&")
                    outfile.write("\n")
                    for data in out_data:
                        outfile.write(data)
                    data_lst = list()
                    product_lst = list()
                    trans_lst = list()
                    protein_lst = list()
                    ID_dic = dict()
                    GENE(Ln[8])
            elif Ln2.startswith("ID=rna"):
                RNA(Ln[8])
            elif Ln2.startswith("ID=cds"):
                CDS(Ln[8])
            else:
                continue
    if data_lst:
        if len(product_lst) > 1:
            Product = ("; ").join(product_lst)
        elif not product_lst:
            Product = "."
        else:
            Product = product_lst[0]
        if len(trans_lst) > 1:
            Trans = ("; ").join(trans_lst)
        elif not trans_lst:
            Trans = "."
        else:
            Trans = trans_lst[0]
        if len(protein_lst) > 1:
            Protein = ("; ").join(protein_lst)
        elif not protein_lst:
            Protein = "."
        else:
            Protein = protein_lst[0]
        data_lst += Product, Trans, Protein
        out_data = (("&\t&").join(data_lst)).strip().split("&")
        outfile.write("\n")
        for data in out_data:
            outfile.write(data)
    infile.close()
    outfile.close()


def GENE(file):
    global data_lst, ID_dic, product_lst, trans_lst, protein_lst
    Ln = file.strip().split(";")
    if [i for i in range(len(Ln)) if "ID=" in Ln[i]]:
        ID = Ln[[i for i in range(len(Ln)) if "ID=" in Ln[i]][0]][
            Ln[[i for i in range(len(Ln)) if "ID=" in Ln[i]][0]].find("=") + 1 :
        ]
        ID_dic["ID"] = ID
    else:
        ID = "."
        ID_dic["ID"] = ID
    if [i for i in range(len(Ln)) if "GeneID" in Ln[i]]:
        GeneID = Ln[[i for i in range(len(Ln)) if "GeneID" in Ln[i]][0]][
            Ln[[i for i in range(len(Ln)) if "GeneID" in Ln[i]][0]].find(
                "GeneID:"
            )
            + 7 :
        ]
    else:
        GeneID = "."
    if [i for i in range(len(Ln)) if "CGNC" in Ln[i]]:
        CGNC = Ln[[i for i in range(len(Ln)) if "CGNC" in Ln[i]][0]][
            Ln[[i for i in range(len(Ln)) if "CGNC" in Ln[i]][0]].find("CGNC:")
            + 5 : Ln[[i for i in range(len(Ln)) if "CGNC" in Ln[i]][0]].find(
                ","
            )
        ]
    else:
        CGNC = "."
    if [i for i in range(len(Ln)) if "gbkey" in Ln[i]]:
        gbkey = Ln[[i for i in range(len(Ln)) if "gbkey" in Ln[i]][0]][
            Ln[[i for i in range(len(Ln)) if "gbkey" in Ln[i]][0]].find("=")
            + 1 :
        ]
    else:
        gbkey = "."
    if [i for i in range(len(Ln)) if "gene_biotype" in Ln[i]]:
        gene_biotype = Ln[
            [i for i in range(len(Ln)) if "gene_biotype" in Ln[i]][0]
        ][
            Ln[[i for i in range(len(Ln)) if "gene_biotype" in Ln[i]][0]].find(
                "="
            )
            + 1 :
        ]
    else:
        gene_biotype = "."
    if [i for i in range(len(Ln)) if "gene" in Ln[i]]:
        gene = Ln[[i for i in range(len(Ln)) if "gene" in Ln[i]][0]][
            Ln[[i for i in range(len(Ln)) if "gene" in Ln[i]][0]].find("=")
            + 1 :
        ]
    else:
        gene = "."
    if [i for i in range(len(Ln)) if "gene_synonym" in Ln[i]]:
        gene_synonym = Ln[
            [i for i in range(len(Ln)) if "gene_synonym" in Ln[i]][0]
        ][
            Ln[[i for i in range(len(Ln)) if "gene_synonym" in Ln[i]][0]].find(
                "="
            )
            + 1 :
        ]
    else:
        gene_synonym = "."
    data_lst += ID, GeneID, CGNC, gbkey, gene_biotype, gene, gene_synonym


def RNA(file):
    global data_lst, ID_dic, product_lst, trans_lst
    Ln = file.strip().split(";")
    Parent = Ln[[i for i in range(len(Ln)) if "Parent" in Ln[i]][0]][
        Ln[[i for i in range(len(Ln)) if "Parent" in Ln[i]][0]].find("=") + 1 :
    ]
    if ID_dic["ID"] == Parent:
        rna_id = Ln[[i for i in range(len(Ln)) if "ID=" in Ln[i]][0]][
            Ln[[i for i in range(len(Ln)) if "ID=" in Ln[i]][0]].find("=") + 1 :
        ]
        ID_dic["RNA"] = rna_id
        if [i for i in range(len(Ln)) if "product" in Ln[i]]:
            product = Ln[[i for i in range(len(Ln)) if "product" in Ln[i]][0]][
                Ln[[i for i in range(len(Ln)) if "product" in Ln[i]][0]].find(
                    "="
                )
                + 1 :
            ]
        else:
            product = "."
        if [i for i in range(len(Ln)) if "transcript_id" in Ln[i]]:
            trans_id = Ln[
                [i for i in range(len(Ln)) if "transcript_id" in Ln[i]][0]
            ][
                Ln[
                    [i for i in range(len(Ln)) if "transcript_id" in Ln[i]][0]
                ].find("=")
                + 1 :
            ]
        else:
            trans_id = "."
        if "%2C" in product:
            product = product.strip().split("%2C")
            for content in product:
                if " X" in content:
                    content = content[: content.find(" X")]
                if content not in product_lst:
                    product_lst.append(content)
        else:
            if " X" in product:
                product = product[: product.find(" X")]
            if product not in product_lst:
                product_lst.append(product)
        if "%2C" in trans_id:
            trans_id = trans_id.strip().split("%2C")
            for content in trans_id:
                if " X" in content:
                    content = content[: content.find(" X")]
                if content not in trans_lst:
                    trans_lst.append(content)
        else:
            if " X" in trans_id:
                trans_id = trans_id[: trans_id.find(" X")]
            if trans_id not in trans_lst:
                trans_lst.append(trans_id)


def CDS(file):
    global data_lst, ID_dic, product_lst, protein_lst
    Ln = file.strip().split(";")
    Parent = Ln[[i for i in range(len(Ln)) if "Parent" in Ln[i]][0]][
        Ln[[i for i in range(len(Ln)) if "Parent" in Ln[i]][0]].find("=") + 1 :
    ]
    try:
        if ID_dic["RNA"] == Parent:
            if [i for i in range(len(Ln)) if "product" in Ln[i]]:
                product = Ln[
                    [i for i in range(len(Ln)) if "product" in Ln[i]][0]
                ][
                    Ln[
                        [i for i in range(len(Ln)) if "product" in Ln[i]][0]
                    ].find("=")
                    + 1 :
                ]
            else:
                product = "."
            if [i for i in range(len(Ln)) if "protein_id" in Ln[i]]:
                protein_id = Ln[
                    [i for i in range(len(Ln)) if "protein_id" in Ln[i]][0]
                ][
                    Ln[
                        [i for i in range(len(Ln)) if "protein_id" in Ln[i]][0]
                    ].find("=")
                    + 1 :
                ]
            else:
                protein_id = "."
            if "%2C" in product:
                product = product.strip().split("%2C")
                for content in product:
                    if " X" in content:
                        content = content[: content.find(" X")]
                    if content not in product_lst:
                        product_lst.append(content)
            else:
                if " X" in product:
                    product = product[: product.find(" X")]
                if product not in product_lst:
                    product_lst.append(product)
            if "%2C" in protein_id:
                protein_id = protein_id.strip().split("%2C")
                for content in protein_id:
                    if " X" in content:
                        content = content[: content.find(" X")]
                    if content not in protein_lst:
                        protein_lst.append(content)
            else:
                if " X" in protein_id:
                    protein_id = protein_id[: protein_id.find(" X")]
                if protein_id not in protein_lst:
                    protein_lst.append(protein_id)
    except:  # V_gene_Segment Case
        if [i for i in range(len(Ln)) if "product" in Ln[i]]:
            product = Ln[[i for i in range(len(Ln)) if "product" in Ln[i]][0]][
                Ln[[i for i in range(len(Ln)) if "product" in Ln[i]][0]].find(
                    "="
                )
                + 1 :
            ]
        else:
            product = "."
        if [i for i in range(len(Ln)) if "protein_id" in Ln[i]]:
            protein_id = Ln[
                [i for i in range(len(Ln)) if "protein_id" in Ln[i]][0]
            ][
                Ln[
                    [i for i in range(len(Ln)) if "protein_id" in Ln[i]][0]
                ].find("=")
                + 1 :
            ]
        else:
            protein_id = "."
        if "%2C" in product:
            product = product.strip().split("%2C")
            for content in product:
                if " X" in content:
                    content = content[: content.find(" X")]
                if content not in product_lst:
                    product_lst.append(content)
        else:
            if " X" in product:
                product = product[: product.find(" X")]
            if product not in product_lst:
                product_lst.append(product)
        if "%2C" in protein_id:
            protein_id = protein_id.strip().split("%2C")
            for content in protein_id:
                if " X" in content:
                    content = content[: content.find(" X")]
                if content not in protein_lst:
                    protein_lst.append(content)
        else:
            if " X" in protein_id:
                protein_id = protein_id[: protein_id.find(" X")]
            if protein_id not in protein_lst:
                protein_lst.append(protein_id)


def main():
    Parsing(inp)


if __name__ == "__main__":
    main()
