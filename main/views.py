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

