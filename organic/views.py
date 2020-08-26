from django.shortcuts import render
from django import forms
from django.http import HttpResponseRedirect

# Defines the organic identification function
def organic_identification(structuralFormula):

    # Defines the formula formatting function which takes in initial structural formula and formats it (case, Cl, etc.)
    def formula_formatting(structuralFormula):

        # Formats the case of the structural formula to be consistent with standard structural formulae
        unformattedStructuralFormula = structuralFormula
        structuralFormula = ""
        for currentChar in unformattedStructuralFormula:
            if currentChar in ["r", "R", "l", "L"]:
                structuralFormula += currentChar.lower()
            else:
                structuralFormula += currentChar.upper()

        # Changes any chlorine atoms (Cl) to 'Xl' to avoid confusion with the 'C' of carbon atoms
        for i in range(1, len(structuralFormula)):
            if structuralFormula[i - 1: i + 1] == "Cl":
                structuralFormula = structuralFormula[: i - 1] + "X" + structuralFormula[i:]

        # Returns an error if there is not a C atom in the structural formula
        if "C" not in structuralFormula:
            return ""
        else:
            return structuralFormula

    # -----------------------------------------------------------------------------------------------------------------

    # Defines the structural separation function that takes in formatted structural formula and splits it into an array
    def structural_separation(structuralFormula):

        # Counts the number of carbon atoms in the structural formula
        carbonAtomPosArray = []
        counterInBracket = False
        for i in range(len(structuralFormula)):
            if structuralFormula[i] == "C" and not counterInBracket:
                carbonAtomPosArray.append(i)

            elif structuralFormula[i] == "(":
                counterInBracket = True

            elif structuralFormula[i] == ")":
                counterInBracket = False

        # Initialises array to hold molecule connections to each main chain carbon atom
        structuralSeparationArray = [["C"] for i in range(len(carbonAtomPosArray))]

        # Loops for each main chain carbon atom and adds corresponding molecules to it's particular array
        for i in range(len(carbonAtomPosArray)):
            if i < len(carbonAtomPosArray) - 1:
                currentMoleculesString = structuralFormula[carbonAtomPosArray[i] + 1: carbonAtomPosArray[i + 1]]
            else:
                currentMoleculesString = structuralFormula[carbonAtomPosArray[i] + 1:]

            # Deals with ketones
            if currentMoleculesString == "O":
                structuralSeparationArray[i].append("O")

            # Deals with OOH, HO, N
            elif currentMoleculesString in ["OOH", "HO", "N"]:
                structuralSeparationArray[i].append(currentMoleculesString)

            # Deals with other other types of molecule(s)
            else:
                counterInBracket = False
                specialGroup = False
                currentSeparatedMolecule = ""
                for j in range(len(currentMoleculesString)):

                    # Deals with non-bracketed molecules
                    if currentMoleculesString[j] in ["H", "F", "I",
                                                     "B", "X"] and not counterInBracket and not specialGroup:
                        if currentMoleculesString[j] in ["H", "F", "I"]:
                            currentSeparatedMolecule = currentMoleculesString[j]
                        else:
                            currentSeparatedMolecule = currentMoleculesString[j: j + 2]
                        structuralSeparationArray[i].append(currentSeparatedMolecule)

                    # Deals with OH/alcohol groups (only at the end of entire structural formula)
                    elif currentMoleculesString[j] == "O" and not counterInBracket:
                        if currentMoleculesString[j + 1] == "H":
                            currentSeparatedMolecule = currentMoleculesString[j: j + 2]
                            specialGroup = True
                            structuralSeparationArray[i].append(currentSeparatedMolecule)

                    # Deals with NH2/amine groups (only at the end of entire structural formula)
                    elif currentMoleculesString[j] == "N" and not counterInBracket:
                        if currentMoleculesString[j + 1: j + 3] == "H2":
                            currentSeparatedMolecule = currentMoleculesString[j: j + 3]
                            specialGroup = True
                            structuralSeparationArray[i].append(currentSeparatedMolecule)

                    # Deals with bracketed molecules
                    elif currentMoleculesString[j] == "(":

                        # Finds current bracketed molecule (in case of multiple brackets in currentMoleculeString)
                        for k in range((len(currentMoleculesString) - 1), j, -1):
                            if currentMoleculesString[k] == ")":
                                closingBracketPos = k

                        currentSeparatedMolecule = currentMoleculesString[j + 1: closingBracketPos]
                        counterInBracket = True
                        structuralSeparationArray[i].append(currentSeparatedMolecule)

                    elif currentMoleculesString[j] == ")":
                        counterInBracket = False

                    # Deals with R2, R3, R4, (R)2, (R)3, (R)4, etc.
                    elif currentMoleculesString[j].isnumeric() and not counterInBracket and not specialGroup:
                        for k in range(int(currentMoleculesString[j]) - 1):
                            structuralSeparationArray[i].append(currentSeparatedMolecule)

        # Deals with formulae with a preceding non-carbon atom molecule
        if structuralFormula[0] != "C":

            currentMoleculesString = structuralFormula[0: carbonAtomPosArray[0]]

            # Deals with non-bracketed molecules
            if len(currentMoleculesString) == 1:
                structuralSeparationArray[0].append(currentMoleculesString)

        # TESTING
        # print("Structural Separation Array", structuralSeparationArray)

        return structuralSeparationArray

    # -----------------------------------------------------------------------------------------------------------------

    # Defines the Atom class
    class Molecule:

        # Defines the class constructor
        def __init__(self, name):
            self.__name = str(name)
            self.__connections = 0
            self.__bonds = []

            # Determines the atom bond capacity based on it's type
            if name == "C":
                self.__capacity = 4
            elif name in ["N", "OOH", "HO"]:
                self.__capacity = 3
            elif name == "O":
                self.__capacity = 2
            elif name in ["H", "F", "I", "Br", "Xl", "OH", "NH2", "CH3"]:
                self.__capacity = 1
            elif len(name) == 4 or len(name) == 5:
                if name[1].isnumeric() and name[3:].isnumeric():
                    if 2 * (int(name[1])) + 1 == int(name[3:]):
                        self.__capacity = 1
            elif len(name) == 6:
                if name[1: 3].isnumeric() and name[4:].isnumeric():
                    if 2 * (int(name[1: 3])) + 1 == int(name[4:]):
                        self.__capacity = 1

        # Defines the get/add methods
        def get_name(self):
            return self.__name

        def get_capacity(self):
            return self.__capacity

        def get_connections(self):
            total = 0
            for currentBond in self.__bonds:
                total += currentBond.get_type()
            return total

        def get_bonds(self):
            return self.__bonds

        def add_bond(self, new):
            self.__bonds.append(new)

    # -----------------------------------------------------------------------------------------------------------------

    # Defines the Bond class
    class Bond:

        # Defines the class constructor
        def __init__(self, molecule1, molecule2):
            self.__molecule1 = molecule1
            self.__molecule2 = molecule2

            # Determines the bond type (if possible) based on the connected molecule names
            if self.__molecule2.get_name() in ["OOH", "N", "HO"]:
                self.__type = 3
            elif self.__molecule2.get_name == "O":
                self.__type = 2
            else:
                self.__type = 1

        # Defines the get/set methods
        def get_type(self):
            return self.__type

        def set_type(self, new):
            self.__type = int(new)

        def get_molecule1(self):
            return self.__molecule1

        def set_molecule1(self, new):
            self.__molecule1 = new

        def get_molecule2(self):
            return self.__molecule2

        def set_molecule2(self, new):
            self.__molecule2 = new

    # -----------------------------------------------------------------------------------------------------------------

    # Defines the structural joining function which converts the separated array into OO format
    def structural_joining(structuralSeparationArray):

        structuralCarbonAtoms = []
        structuralOtherMolecules = []
        structuralBonds = []

        # Creates a Molecule objects for each molecule in the entire structure
        for i in range(len(structuralSeparationArray)):
            structuralCarbonAtoms.append(Molecule(structuralSeparationArray[i][0]))

            for j in range(1, len(structuralSeparationArray[i])):
                structuralOtherMolecules.append(Molecule(structuralSeparationArray[i][j]))

        # Creates bonds between each main chain carbon and their corresponding molecule branches
        currentBondCounter = 0
        for i in range(len(structuralSeparationArray)):
            for j in range(1, len(structuralSeparationArray[i])):

                structuralBonds.append(Bond(structuralCarbonAtoms[i], structuralOtherMolecules[currentBondCounter]))
                structuralCarbonAtoms[i].add_bond(structuralBonds[-1])
                structuralOtherMolecules[currentBondCounter].add_bond(structuralBonds[-1])
                currentBondCounter += 1

        # Creates bonds between the main chain carbon atoms themselves
        for i in range(len(structuralCarbonAtoms) - 1):
            structuralBonds.append(Bond(structuralCarbonAtoms[i], structuralCarbonAtoms[i + 1]))
            structuralCarbonAtoms[i].add_bond(structuralBonds[-1])
            structuralCarbonAtoms[i + 1].add_bond(structuralBonds[-1])

        # Changes the type of carbon-other bonds to find double and triple bonds accordingly (within capacity of atoms)
        bondTypeIncrements = -1
        while bondTypeIncrements != 0:
            bondTypeIncrements = 0
            for i in range(len(structuralBonds)):
                currentMolecule1 = structuralBonds[i].get_molecule1()
                currentMolecule2 = structuralBonds[i].get_molecule2()

                if currentMolecule1.get_capacity() - currentMolecule1.get_connections() > 0 and \
                        currentMolecule2.get_capacity() - currentMolecule2.get_connections() > 0 and \
                        structuralBonds[i].get_type() < 3:

                    structuralBonds[i].set_type(structuralBonds[i].get_type() + 1)
                    bondTypeIncrements += 1

        # TESTING
        # print("C Atoms", [structuralCarbonAtoms[i].get_name() for i in range(len(structuralCarbonAtoms))])
        # print("Molecules", [structuralOtherMolecules[i].get_name() for i in range(len(structuralOtherMolecules))])
        # print("Bond Molecule 1", [structuralBonds[i].get_molecule1().get_name() for i in range(len(structuralBonds))])
        # print("Bond Molecule 2", [structuralBonds[i].get_molecule2().get_name() for i in range(len(structuralBonds))])
        # print("Bond Type", [structuralBonds[i].get_type() for i in range(len(structuralBonds))])

        return [structuralCarbonAtoms, structuralOtherMolecules, structuralBonds]

    # -----------------------------------------------------------------------------------------------------------------

    # Defines the structural validation function which ensures each atom has reached their exact bond capacity
    def structural_validation(structuralJoiningArray):
        structureValid = True

        # Validates the carbon atoms
        try:
            for currentAtom in structuralJoiningArray[0]:
                if currentAtom.get_capacity() != currentAtom.get_connections():

                    structureValid = False

        except AttributeError:
            structureValid = False

        # Validates the other molecules
        try:
            for currentMolecule in structuralJoiningArray[1]:
                if currentMolecule.get_capacity() != currentMolecule.get_connections():

                    structureValid = False

        except AttributeError:
            structureValid = False

        return structureValid

    # -----------------------------------------------------------------------------------------------------------------

    # Defines the organic name creation function which takes in the OO format of the structural compound
    def organic_name_creation(structuralJoiningArray):

        structuralCarbonAtoms = structuralJoiningArray[0]

        # Initialises arrays with IUPAC nomenclature names/prefixes
        parentChainNamesArray = ["", "meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec",
                                 "undec", "dodec", "tridec", "tetradec", "pentadec", "hexadec", "heptadec", "octadec",
                                 "nonadec", "icos"]
        multipleChainNamesArray = ["", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca",
                                   "undeca", "dodeca", "trideca", "tetradeca", "pentadeca", "hexadeca", "heptadeca",
                                   "octadeca", "nonadeca", "icosa"]
        groupPrecedenceArray = ["OOH", "N", "HO", "O", "OH", "NH2"]
        groupPrefixArray = ["carboxy", "cyano", "oxo", "oxo", "hydroxy", "amino"]
        groupSuffixArray = ["oic acid", "nitrile", "al", "one", "ol", "amine"]
        otherMoleculeNamesArray = ["Br", "Xl", "F", "I", "CH3"]
        otherMoleculesPrefixArray = ["bromo", "chloro", "fluoro", "iodo", "methyl"]

        # Initialises arrays/variables for storing information about the structural compound
        prominentFunctionalGroup = ["initial", 0, []]
        otherFunctionalGroupsArray = []
        otherMoleculeBranchesArray = []
        carbonMultipleBondsArray = [["en", 0, []], ["yn", 0, []]]
        carbonAlkeneAlkyneBondsArray = []

        # Loops through each main chain carbon atom and adds functional groups to an array appropriately
        for i in range(len(structuralCarbonAtoms)):
            for currentBond in structuralCarbonAtoms[i].get_bonds():

                currentMolecule = currentBond.get_molecule2()

                # Deals with functional groups with precedence
                if currentMolecule.get_name() in groupPrecedenceArray:

                    # Deals with the first functional group of this type
                    if prominentFunctionalGroup[0] == "initial":
                        prominentFunctionalGroup = [currentMolecule.get_name(), 1, [i + 1]]

                    # Deals with functional groups of greater precedence
                    elif groupPrecedenceArray.index(prominentFunctionalGroup[0]) > \
                            groupPrecedenceArray.index(currentMolecule.get_name()):
                        otherFunctionalGroupsArray.append(prominentFunctionalGroup)
                        prominentFunctionalGroup = [currentMolecule.get_name(), 1, [i + 1]]

                    # Deals with functional groups of equal precedence
                    elif groupPrecedenceArray.index(prominentFunctionalGroup[0]) == \
                            groupPrecedenceArray.index(currentMolecule.get_name()):
                        prominentFunctionalGroup[1] += 1
                        prominentFunctionalGroup[2].append(i + 1)

                    # Deals with functional groups of lower precedence
                    else:
                        newFunctionalGroup = True

                        for j in range(len(otherFunctionalGroupsArray)):
                            if otherFunctionalGroupsArray[j][0] == currentMolecule.get_name():
                                otherFunctionalGroupsArray[j][1] += 1
                                otherFunctionalGroupsArray[j][2].append(i + 1)
                                newFunctionalGroup = False

                        if newFunctionalGroup:
                            otherFunctionalGroupsArray.append([currentMolecule.get_name(), 1, [i + 1]])

                # Deals with other molecule branches (F, I, Br, Xl)
                elif currentMolecule.get_name() in otherMoleculeNamesArray:
                    newMoleculeBranch = True

                    for j in range(len(otherMoleculeBranchesArray)):
                        if otherMoleculeBranchesArray[j][0] == currentMolecule.get_name():

                            otherMoleculeBranchesArray[j][1] += 1
                            otherMoleculeBranchesArray[j][2].append(i + 1)
                            newMoleculeBranch = False

                    if newMoleculeBranch:
                        otherMoleculeBranchesArray.append([currentMolecule.get_name(), 1, [i + 1]])

                # Deals with other alkyl molecule branches (C2H5, C3H7, C4H9, C5H11, etc.)
                elif 4 <= len(currentMolecule.get_name()) <= 5:
                    if currentMolecule.get_name()[1].isnumeric() and currentMolecule.get_name()[3:].isnumeric():
                        if int(currentMolecule.get_name()[1]) * 2 + 1 == int(currentMolecule.get_name()[3:]):
                            newMoleculeBranch = True

                            for j in range(len(otherMoleculeBranchesArray)):
                                if otherMoleculeBranchesArray[j][0] == currentMolecule.get_name():
                                    otherMoleculeBranchesArray[j][1] += 1
                                    otherMoleculeBranchesArray[j][2].append(i + 1)
                                    newMoleculeBranch = False

                            if newMoleculeBranch:
                                otherMoleculeBranchesArray.append([currentMolecule.get_name(), 1, [i + 1]])

                # Deals with other alkyl molecule branches (C10H21, C11H23, etc.)
                elif len(currentMolecule.get_name()) == 6:
                    if currentMolecule.get_name()[1: 3].isnumeric() and currentMolecule.get_name()[4:].isnumeric():
                        if int(currentMolecule.get_name()[1: 3]) * 2 + 1 == int(currentMolecule.get_name()[4:]):
                            newMoleculeBranch = True

                            for j in range(len(otherMoleculeBranchesArray)):
                                if otherMoleculeBranchesArray[j][0] == currentMolecule.get_name():
                                    otherMoleculeBranchesArray[j][1] += 1
                                    otherMoleculeBranchesArray[j][2].append(i + 1)
                                    newMoleculeBranch = False

                            if newMoleculeBranch:
                                otherMoleculeBranchesArray.append([currentMolecule.get_name(), 1, [i + 1]])

                # Deals with alkyne/alkene indicators (i.e. double and triple carbon bonds)
                elif currentMolecule.get_name() == "C" and currentBond.get_type() in [2, 3] and \
                        currentBond not in carbonAlkeneAlkyneBondsArray:
                    carbonMultipleBondsArray[currentBond.get_type() - 2][1] += 1
                    carbonMultipleBondsArray[currentBond.get_type() - 2][2].append(i + 1)
                    carbonAlkeneAlkyneBondsArray.append(currentBond)

        # TESTING
        # print("Prominent Functional Group", prominentFunctionalGroup)
        # print("Other Functional Groups", otherFunctionalGroupsArray)
        # print("Other Molecule Branches", otherMoleculeBranchesArray)
        # print("Alkene/Alkyne Counter", carbonMultipleBondsArray)

        # Creates an array with the names of the prefixes to add (unsorted)
        prefixesToAddArray = otherFunctionalGroupsArray[:] + otherMoleculeBranchesArray[:]

        for i in range(len(otherFunctionalGroupsArray)):
            prefixesToAddArray[i][0] = groupPrefixArray[groupPrecedenceArray.index(otherFunctionalGroupsArray[i][0])]

        for j in range(len(otherFunctionalGroupsArray), len(otherFunctionalGroupsArray) + len(otherMoleculeBranchesArray)):
            if prefixesToAddArray[j][0] in otherMoleculeNamesArray:
                prefixesToAddArray[j][0] = otherMoleculesPrefixArray[otherMoleculeNamesArray.index(prefixesToAddArray[j][0])]

            elif prefixesToAddArray[j][0][1: 3].isnumeric():
                prefixesToAddArray[j][0] = parentChainNamesArray[int(prefixesToAddArray[j][0][1: 3])] + "yl"
            elif prefixesToAddArray[j][0][1].isnumeric():
                prefixesToAddArray[j][0] = parentChainNamesArray[int(prefixesToAddArray[j][0][1])] + "yl"

        # Sorts the 'prefixes to add' array in alphabetical order
        unsortedPrefixCount = len(prefixesToAddArray)
        swapCount = -1
        while swapCount != 0:
            swapCount = 0

            for i in range(0, unsortedPrefixCount - 1):

                if prefixesToAddArray[i][0] > prefixesToAddArray[i + 1][0]:
                    tempPrefix = prefixesToAddArray[i][:]
                    prefixesToAddArray[i] = prefixesToAddArray[i + 1][:]
                    prefixesToAddArray[i + 1] = tempPrefix[:]
                    swapCount += 1

            unsortedPrefixCount -= 1

        # Decides position numbering based on prominent group and minimising the total of the position values
        carbonChainLength = len(structuralCarbonAtoms)
        normalPositionTotal = 0
        inversePositionTotal = 0
        for i in range(len(prefixesToAddArray)):
            for j in range(prefixesToAddArray[i][1]):
                normalPositionTotal += prefixesToAddArray[i][2][j]
                inversePositionTotal += carbonChainLength + 1 - prefixesToAddArray[i][2][j]

        for i in range(len(carbonMultipleBondsArray)):
            for j in range(carbonMultipleBondsArray[i][1]):
                normalPositionTotal += carbonMultipleBondsArray[i][2][j]
                inversePositionTotal += carbonChainLength - carbonMultipleBondsArray[i][2][j]

        for i in range(prominentFunctionalGroup[1]):
            normalPositionTotal += prominentFunctionalGroup[2][i]
            inversePositionTotal += carbonChainLength + 1 - prominentFunctionalGroup[2][i]

        # Decides on whether to switch the position indexes based on multiple conditions
        switchPositionIndexes = False
        if prominentFunctionalGroup[1] != 0:
            if prominentFunctionalGroup[1] == prominentFunctionalGroup[2][0] == 1:
                switchPositionIndexes = False
            elif prominentFunctionalGroup[1] == 1 and prominentFunctionalGroup[2][0] == carbonChainLength != 1:
                switchPositionIndexes = True
            elif inversePositionTotal < normalPositionTotal:
                switchPositionIndexes = True

        elif inversePositionTotal < normalPositionTotal:
            switchPositionIndexes = True

        # Switches position indexes to start from the opposite end of the carbon chain
        if switchPositionIndexes:
            for i in range(len(prefixesToAddArray)):
                for j in range(prefixesToAddArray[i][1]):
                    prefixesToAddArray[i][2][j] = carbonChainLength + 1 - prefixesToAddArray[i][2][j]

            for i in range(len(carbonMultipleBondsArray)):
                for j in range(carbonMultipleBondsArray[i][1]):

                    carbonMultipleBondsArray[i][2][j] = carbonChainLength - carbonMultipleBondsArray[i][2][j]

            for i in range(prominentFunctionalGroup[1]):
                prominentFunctionalGroup[2][i] = carbonChainLength + 1 - prominentFunctionalGroup[2][i]

        # Sorts the position list for each prefix and double/triple bond counters
        for i in range(len(prefixesToAddArray)):
            prefixesToAddArray[i][2].sort()
        carbonMultipleBondsArray[0][2].sort()
        carbonMultipleBondsArray[1][2].sort()

        # Begins to construct the organic compound name based on the gathered information
        organicCompoundName = ""

        # Adds the prefixes from the 'prefixes to add' array into the organic compound name
        for i in range(len(prefixesToAddArray)):
            for j in range(prefixesToAddArray[i][1]):
                organicCompoundName += str(prefixesToAddArray[i][2][j]) + ","
            organicCompoundName = organicCompoundName[: -1] + "-"
            organicCompoundName += multipleChainNamesArray[prefixesToAddArray[i][1]] + prefixesToAddArray[i][0] + "-"

        # Adds the main chain parent name to the organic compound name
        organicCompoundName += parentChainNamesArray[carbonChainLength]

        # Adds the alkyne/alkene names to the organic compound name
        if carbonMultipleBondsArray[1][1] != 0:
            organicCompoundName += "-"
            for i in range(carbonMultipleBondsArray[1][1]):
                organicCompoundName += str(carbonMultipleBondsArray[1][2][i]) + ","
            organicCompoundName = organicCompoundName[: -1] + "-"
            organicCompoundName += multipleChainNamesArray[carbonMultipleBondsArray[1][1]]
            organicCompoundName += carbonMultipleBondsArray[1][0]

        elif carbonMultipleBondsArray[0][1] != 0:
            organicCompoundName += "-"
            for i in range(carbonMultipleBondsArray[0][1]):
                organicCompoundName += str(carbonMultipleBondsArray[0][2][i]) + ","
            organicCompoundName = organicCompoundName[: -1] + "-"
            organicCompoundName += multipleChainNamesArray[carbonMultipleBondsArray[0][1]]
            organicCompoundName += carbonMultipleBondsArray[0][0]

        else:
            organicCompoundName += "an"

        # Adds the organic compound name ending using the prominent functional group
        if prominentFunctionalGroup[1] != 0:
            organicCompoundName += "-"
            for i in range(prominentFunctionalGroup[1]):
                organicCompoundName += str(prominentFunctionalGroup[2][i]) + ","
            organicCompoundName = organicCompoundName[: -1] + "-" + multipleChainNamesArray[prominentFunctionalGroup[1]]

        else:
            organicCompoundName += multipleChainNamesArray[prominentFunctionalGroup[1]]

        if prominentFunctionalGroup[0] in groupPrecedenceArray:
            organicCompoundName += groupSuffixArray[groupPrecedenceArray.index(prominentFunctionalGroup[0])]

        else:
            organicCompoundName += "e"

        # Edits the organic compound name for specific amine, nitrile and carboxylic acid wordings
        if organicCompoundName[-5:] == "amine":
            organicCompoundName = organicCompoundName.replace("an", "yl")
        if organicCompoundName[-7:] == "nitrile":
            organicCompoundName = organicCompoundName.replace("an", "ane")
            organicCompoundName = organicCompoundName.replace("en", "ene")
            organicCompoundName = organicCompoundName.replace("yn", "yne")

        return organicCompoundName

    # -----------------------------------------------------------------------------------------------------------------

    # Completes the organic identification by combining each of the subroutines
    formattedStructuralFormula = formula_formatting(structuralFormula)

    if formattedStructuralFormula == "":
        return "unrecognisable by the program. Please ensure there is at least one C atom."

    else:
        structuralSeparationArray = structural_separation(formattedStructuralFormula)
        structuralJoiningArray = structural_joining(structuralSeparationArray)

        if structural_validation(structuralJoiningArray):
            return organic_name_creation(structuralJoiningArray)
        else:
            return "unrecognisable by the program. Please try another formula."


# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------


# Creates an organic compound input form
class OrganicForm(forms.Form):
    structuralFormula = forms.CharField(label="Structural Formula")


# Create your views here.
def index(request):
    form = OrganicForm()
    return render(request, "organic/index.html", {"form": form})


def organic_identifier(request, chemicalFormula):

    chemicalName = organic_identification(chemicalFormula)
    return render(request, "organic/identification.html", {
        "chemicalName": chemicalName,
        "chemicalFormula": chemicalFormula})


def organic_identifier_form(request):

    if request.method == "POST":
        form = OrganicForm(request.POST)

        if form.is_valid():
            chemicalFormula = form.cleaned_data["structuralFormula"]
            return HttpResponseRedirect("/organic/formula/" + chemicalFormula)

    else:
        form = OrganicForm()

    return render(request, "index.html", {"form": form})