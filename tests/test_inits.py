
import pytest
from gastronomy.molecules import Molecule
from gastronomy.minerals import Mineral

def test_molecules():
    result = Molecule({'C':1, 'O':1})
    assert isinstance(result, Molecule)

def test_minerals():
    result = Mineral({'Si':1, 'O':2})
    assert isinstance(result, Mineral)


