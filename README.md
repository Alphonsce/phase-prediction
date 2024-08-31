# Качественное предсказание фазового состояния органических веществ при нормальных условиях. (Жидкость или нет данное вещество)

## Данные:

- Разделитель - пробелы.
- Формат: номер строки, SMILES, CAS, метка жидкость/твердое тело (1 - твердое тело, 2 - жидкость)

CAS - это уникальный идентификатор вещества (без физического смысла). Если упрощенно это его название в виде чисел. Если его вбить в поиск например на https://pubchem.ncbi.nlm.nih.gov/ то оно находит это вещество

## rdkit lib

работа с библиотекой rdkit выглядит примерно так

```
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles(smiles_line) # создаем обьект молекула по SMILES формуле (2 строка датасета)
h_bond_donors= Chem.rdMolDescriptors.CalcNumHBD(m) # считаем число доноров водородной связи в этой молекуле 
```

В этой библиотеке довольно много разных молекулярных дескрипторов https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html

для начала имеет смысл все это запустить и посмотреть например как зависит классификация от количества тяжелых (те всех кроме атомов водорода) атомов в молекуле 

```
 rdkit.Chem.rdMolDescriptors.CalcNumHeavyAtoms((Mol)mol) → int 
    returns the number of heavy atoms for a molecule
```