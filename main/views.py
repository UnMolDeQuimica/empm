from django.shortcuts import render
from django.views.generic import TemplateView, View
from .models import MoleculeRecord
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from django.http import HttpResponseServerError
from rdkit.Chem import Draw, AllChem, rdFMCS, DataStructs, MACCSkeys
from rdkit.Chem.Draw import rdMolDraw2D


class Home(TemplateView):
    template_name = "index.html"


class HowDoesItWork(TemplateView):
    template_name = "how_does_it_work.html"


class CookiesPolicy(TemplateView):
    template_name = "cookies_policy.html"


class Molecule(View):
    template_name = "single_molecule_base.html"
    post_template_name = "molecule.html"
    incorrect_smiles_template = "incorrect_smiles_message.html"
    error_template = "unexpected_error.html"
    empty_template = "empty_smiles_message.html"

    def get(self, request, *args, **kwargs):
        return render(request=request, template_name=self.template_name)

    def post(self, request, *args, **kwargs):
        try:
            smiles = self.request.POST.get("smilesInput")

            if not smiles:
                return render(request=request, template_name=self.empty_template)

            is_valid_smiles = Chem.MolFromSmiles(smiles, sanitize=False)
            if not is_valid_smiles:
                return render(request=request, template_name=self.error_template)

            molecule = self.get_object(smiles=smiles)
            context_data = self.get_context_data(molecule=molecule)

            return render(
                request=request,
                template_name=self.post_template_name,
                context=context_data,
            )

        except Exception as e:
            print(e)
            return HttpResponseServerError()

    def get_object(self, smiles):
        return MoleculeRecord.objects.get_or_create(smiles=smiles)[0]

    def get_context_data(self, molecule: MoleculeRecord):
        context_data = {"molecule": molecule}
        mol = Chem.MolFromSmiles(molecule.smiles)
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().clearBackground = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        context_data["molecule_image"] = drawer.GetDrawingText()
        context_data["mol_block"] = Chem.MolToMolBlock(mol)
        context_data["smiles"] = molecule.smiles
        context_data["maccs"] = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
        context_data["canonical_smiles"] = Chem.MolToCXSmiles(mol)
        context_data["inchi"] = Chem.MolToInchi(mol)
        context_data["inchi_key"] = Chem.MolToInchiKey(mol)
        context_data["smarts"] = Chem.MolToSmarts(mol)
        context_data["molecular_formula"] = rdMolDescriptors.CalcMolFormula(mol)
        context_data["molecular_weight"] = rdMolDescriptors.CalcExactMolWt(mol)
        return context_data


class MoleculeComparison(View):
    template_name = "molecule_comparison_base.html"
    post_template_name = "molecule_comparison.html"
    incorrect_smiles_template = "incorrect_smiles_message.html"
    error_template = "unexpected_error.html"
    empty_template = "empty_smiles_message.html"

    def get(self, request, *args, **kwargs):
        return render(request=request, template_name=self.template_name)

    def post(self, request, *args, **kwargs):
        try:
            smiles_a = self.request.POST.get("smilesInputA")
            smiles_b = self.request.POST.get("smilesInputB")

            if (not smiles_a) or (not smiles_b):
                return render(request=request, template_name=self.empty_template)

            is_valid_smiles = Chem.MolFromSmiles(
                smiles_a, sanitize=False
            ) and Chem.MolFromSmiles(smiles_b, sanitize=False)
            if not is_valid_smiles:
                return render(request=request, template_name=self.error_template)

            molecule_a = self.get_object(smiles=smiles_a)
            molecule_b = self.get_object(smiles=smiles_b)
            context_data = self.get_context_data(
                molecule_a=molecule_a, molecule_b=molecule_b
            )

            return render(
                request=request,
                template_name=self.post_template_name,
                context=context_data,
            )

        except Exception as e:
            print(e)
            return HttpResponseServerError()

    def get_object(self, smiles):
        return MoleculeRecord.objects.get_or_create(smiles=smiles)[0]

    def get_context_data(self, molecule_a: MoleculeRecord, molecule_b: MoleculeRecord):
        context_data = {}

        mol_a = Chem.MolFromSmiles(molecule_a.smiles)
        mol_b = Chem.MolFromSmiles(molecule_b.smiles)
        mcs = rdFMCS.FindMCS([mol_a, mol_b])
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

        context_data["molecule_a"] = molecule_a
        context_data["molecule_b"] = molecule_b
        context_data["molecule_a_svg"] = self.draw_molecule_svg(mol_a, mcs_mol)
        context_data["molecule_b_svg"] = self.draw_molecule_svg(mol_b, mcs_mol)
        context_data = self.calculate_similarities(mol_a, mol_b, context_data)

        return context_data

    def draw_molecule_svg(self, mol: Chem.Mol, mcs_mol: Chem.Mol):
        matches = sum(mol.GetSubstructMatches(mcs_mol), ())
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().clearBackground = False
        drawer.DrawMolecule(mol, highlightAtoms=matches)
        drawer.FinishDrawing()

        return drawer.GetDrawingText()

    def calculate_similarities(self, mol_a, mol_b, context_data):
        morgan1 = AllChem.GetMorganFingerprintAsBitVect(mol_a, 2, nBits=2048)
        morgan2 = AllChem.GetMorganFingerprintAsBitVect(mol_b, 2, nBits=2048)
        tanimoto_morgan = DataStructs.TanimotoSimilarity(morgan1, morgan2)

        fp1 = Chem.RDKFingerprint(mol_a)
        fp2 = Chem.RDKFingerprint(mol_b)
        tanimoto_rdkit = DataStructs.TanimotoSimilarity(fp1, fp2)

        maccs1 = MACCSkeys.GenMACCSKeys(mol_a)
        maccs2 = MACCSkeys.GenMACCSKeys(mol_b)
        tanimoto_maccs = DataStructs.TanimotoSimilarity(maccs1, maccs2)

        dice_sim = DataStructs.DiceSimilarity(morgan1, morgan2)

        cosine_sim = DataStructs.CosineSimilarity(morgan1, morgan2)

        context_data["tanimoto_morgan"] = round(tanimoto_morgan, 4) * 100
        context_data["tanimoto_rdkit"] = round(tanimoto_rdkit, 4) * 100
        context_data["tanimoto_maccs"] = round(tanimoto_maccs, 4) * 100
        context_data["dice"] = round(dice_sim, 4) * 100
        context_data["cosine"] = round(cosine_sim, 4) * 100

        return context_data
