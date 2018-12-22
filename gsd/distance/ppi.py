from pandas import read_table, DataFrame


def load_ppi_mitab(ppi_file: str) -> DataFrame:
    ppi_data = read_table(ppi_file)
    ppi_data = ppi_data.assign(
        FromId=[int(e.split(":")[1]) for e in ppi_data["#ID Interactor A"]],
        ToId=[int(e.split(":")[1]) for e in ppi_data["ID Interactor B"]],
        FromTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor A"]],
        ToTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor B"]])
    return ppi_data
