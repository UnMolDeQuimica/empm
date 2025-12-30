from django.db import models
from admet_ai import ADMETModel


# Create your models here.
class MoleculeRecord(models.Model):
    smiles = models.TextField(verbose_name="SMILES", primary_key=True)
    molecular_weight = models.FloatField(
        verbose_name="molecular_weight", null=True, blank=True
    )
    logp = models.FloatField(verbose_name="logP", null=True, blank=True)
    hydrogen_bond_acceptors = models.FloatField(
        verbose_name="hydrogen_bond_acceptors", null=True, blank=True
    )
    hydrogen_bond_donors = models.FloatField(
        verbose_name="hydrogen_bond_donors", null=True, blank=True
    )
    lipinski = models.FloatField(verbose_name="Lipinski", null=True, blank=True)
    qed = models.FloatField(verbose_name="QED", null=True, blank=True)
    stereo_centers = models.FloatField(
        verbose_name="stereo_centers", null=True, blank=True
    )
    tpsa = models.FloatField(verbose_name="tpsa", null=True, blank=True)
    ames = models.FloatField(verbose_name="AMES", null=True, blank=True)
    bbb_martins = models.FloatField(verbose_name="BBB_Martins", null=True, blank=True)
    bioavailability_ma = models.FloatField(
        verbose_name="Bioavailability_Ma", null=True, blank=True
    )
    cyp1a2_veith = models.FloatField(verbose_name="CYP1A2_Veith", null=True, blank=True)
    cyp2c19_veith = models.FloatField(
        verbose_name="CYP2C19_Veith", null=True, blank=True
    )
    cyp2c9_substrate_carbonmangels = models.FloatField(
        verbose_name="CYP2C9_Substrate_CarbonMangels", null=True, blank=True
    )
    cyp2c9_veith = models.FloatField(verbose_name="CYP2C9_Veith", null=True, blank=True)
    cyp2d6_substrate_carbonmangels = models.FloatField(
        verbose_name="CYP2D6_Substrate_CarbonMangels", null=True, blank=True
    )
    cyp2d6_veith = models.FloatField(verbose_name="CYP2D6_Veith", null=True, blank=True)
    cyp3a4_substrate_carbonmangels = models.FloatField(
        verbose_name="CYP3A4_Substrate_CarbonMangels", null=True, blank=True
    )
    cyp3a4_veith = models.FloatField(verbose_name="CYP3A4_Veith", null=True, blank=True)
    carcinogens_lagunin = models.FloatField(
        verbose_name="Carcinogens_Lagunin", null=True, blank=True
    )
    clintox = models.FloatField(verbose_name="ClinTox", null=True, blank=True)
    dili = models.FloatField(verbose_name="DILI", null=True, blank=True)
    hia_hou = models.FloatField(verbose_name="HIA_Hou", null=True, blank=True)
    nr_ar_lbd = models.FloatField(verbose_name="NR-AR-LBD", null=True, blank=True)
    nr_ar = models.FloatField(verbose_name="NR-AR", null=True, blank=True)
    nr_ahr = models.FloatField(verbose_name="NR-AhR", null=True, blank=True)
    nr_aromatase = models.FloatField(verbose_name="NR-Aromatase", null=True, blank=True)
    nr_er_lbd = models.FloatField(verbose_name="NR-ER-LBD", null=True, blank=True)
    nr_er = models.FloatField(verbose_name="NR-ER", null=True, blank=True)
    nr_ppar_gamma = models.FloatField(
        verbose_name="NR-PPAR-gamma", null=True, blank=True
    )
    pampa_ncats = models.FloatField(verbose_name="PAMPA_NCATS", null=True, blank=True)
    pgp_broccatelli = models.FloatField(
        verbose_name="Pgp_Broccatelli", null=True, blank=True
    )
    sr_are = models.FloatField(verbose_name="SR-ARE", null=True, blank=True)
    sr_atad5 = models.FloatField(verbose_name="SR-ATAD5", null=True, blank=True)
    sr_hse = models.FloatField(verbose_name="SR-HSE", null=True, blank=True)
    sr_mmp = models.FloatField(verbose_name="SR-MMP", null=True, blank=True)
    sr_p53 = models.FloatField(verbose_name="SR-p53", null=True, blank=True)
    skin_reaction = models.FloatField(
        verbose_name="Skin_Reaction", null=True, blank=True
    )
    herg = models.FloatField(verbose_name="hERG", null=True, blank=True)
    caco2_wang = models.FloatField(verbose_name="Caco2_Wang", null=True, blank=True)
    clearance_hepatocyte_az = models.FloatField(
        verbose_name="Clearance_Hepatocyte_AZ", null=True, blank=True
    )
    clearance_microsome_az = models.FloatField(
        verbose_name="Clearance_Microsome_AZ", null=True, blank=True
    )
    half_life_obach = models.FloatField(
        verbose_name="Half_Life_Obach", null=True, blank=True
    )
    hydrationfreeenergy_freesolv = models.FloatField(
        verbose_name="HydrationFreeEnergy_FreeSolv", null=True, blank=True
    )
    ld50_zhu = models.FloatField(verbose_name="LD50_Zhu", null=True, blank=True)
    lipophilicity_astrazeneca = models.FloatField(
        verbose_name="Lipophilicity_AstraZeneca", null=True, blank=True
    )
    ppbr_az = models.FloatField(verbose_name="PPBR_AZ", null=True, blank=True)
    solubility_aqsoldb = models.FloatField(
        verbose_name="Solubility_AqSolDB", null=True, blank=True
    )
    vdss_lombardo = models.FloatField(
        verbose_name="VDss_Lombardo", null=True, blank=True
    )
    molecular_weight_drugbank_approved_percentile = models.FloatField(
        verbose_name="molecular_weight_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    logp_drugbank_approved_percentile = models.FloatField(
        verbose_name="logP_drugbank_approved_percentile", null=True, blank=True
    )
    hydrogen_bond_acceptors_drugbank_approved_percentile = models.FloatField(
        verbose_name="hydrogen_bond_acceptors_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    hydrogen_bond_donors_drugbank_approved_percentile = models.FloatField(
        verbose_name="hydrogen_bond_donors_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    lipinski_drugbank_approved_percentile = models.FloatField(
        verbose_name="Lipinski_drugbank_approved_percentile", null=True, blank=True
    )
    qed_drugbank_approved_percentile = models.FloatField(
        verbose_name="QED_drugbank_approved_percentile", null=True, blank=True
    )
    stereo_centers_drugbank_approved_percentile = models.FloatField(
        verbose_name="stereo_centers_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    tpsa_drugbank_approved_percentile = models.FloatField(
        verbose_name="tpsa_drugbank_approved_percentile", null=True, blank=True
    )
    ames_drugbank_approved_percentile = models.FloatField(
        verbose_name="AMES_drugbank_approved_percentile", null=True, blank=True
    )
    bbb_martins_drugbank_approved_percentile = models.FloatField(
        verbose_name="BBB_Martins_drugbank_approved_percentile", null=True, blank=True
    )
    bioavailability_ma_drugbank_approved_percentile = models.FloatField(
        verbose_name="Bioavailability_Ma_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    cyp1a2_veith_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP1A2_Veith_drugbank_approved_percentile", null=True, blank=True
    )
    cyp2c19_veith_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP2C19_Veith_drugbank_approved_percentile", null=True, blank=True
    )
    cyp2c9_substrate_carbonmangels_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP2C9_Substrate_CarbonMangels_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    cyp2c9_veith_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP2C9_Veith_drugbank_approved_percentile", null=True, blank=True
    )
    cyp2d6_substrate_carbonmangels_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP2D6_Substrate_CarbonMangels_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    cyp2d6_veith_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP2D6_Veith_drugbank_approved_percentile", null=True, blank=True
    )
    cyp3a4_substrate_carbonmangels_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP3A4_Substrate_CarbonMangels_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    cyp3a4_veith_drugbank_approved_percentile = models.FloatField(
        verbose_name="CYP3A4_Veith_drugbank_approved_percentile", null=True, blank=True
    )
    carcinogens_lagunin_drugbank_approved_percentile = models.FloatField(
        verbose_name="Carcinogens_Lagunin_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    clintox_drugbank_approved_percentile = models.FloatField(
        verbose_name="ClinTox_drugbank_approved_percentile", null=True, blank=True
    )
    dili_drugbank_approved_percentile = models.FloatField(
        verbose_name="DILI_drugbank_approved_percentile", null=True, blank=True
    )
    hia_hou_drugbank_approved_percentile = models.FloatField(
        verbose_name="HIA_Hou_drugbank_approved_percentile", null=True, blank=True
    )
    nr_ar_lbd_drugbank_approved_percentile = models.FloatField(
        verbose_name="NR-AR-LBD_drugbank_approved_percentile", null=True, blank=True
    )
    nr_ar_drugbank_approved_percentile = models.FloatField(
        verbose_name="NR-AR_drugbank_approved_percentile", null=True, blank=True
    )
    nr_ahr_drugbank_approved_percentile = models.FloatField(
        verbose_name="NR-AhR_drugbank_approved_percentile", null=True, blank=True
    )
    nr_aromatase_drugbank_approved_percentile = models.FloatField(
        verbose_name="NR-Aromatase_drugbank_approved_percentile", null=True, blank=True
    )
    nr_er_lbd_drugbank_approved_percentile = models.FloatField(
        verbose_name="NR-ER-LBD_drugbank_approved_percentile", null=True, blank=True
    )
    nr_er_drugbank_approved_percentile = models.FloatField(
        verbose_name="NR-ER_drugbank_approved_percentile", null=True, blank=True
    )
    nr_ppar_gamma_drugbank_approved_percentile = models.FloatField(
        verbose_name="NR-PPAR-gamma_drugbank_approved_percentile", null=True, blank=True
    )
    pampa_ncats_drugbank_approved_percentile = models.FloatField(
        verbose_name="PAMPA_NCATS_drugbank_approved_percentile", null=True, blank=True
    )
    pgp_broccatelli_drugbank_approved_percentile = models.FloatField(
        verbose_name="Pgp_Broccatelli_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    sr_are_drugbank_approved_percentile = models.FloatField(
        verbose_name="SR-ARE_drugbank_approved_percentile", null=True, blank=True
    )
    sr_atad5_drugbank_approved_percentile = models.FloatField(
        verbose_name="SR-ATAD5_drugbank_approved_percentile", null=True, blank=True
    )
    sr_hse_drugbank_approved_percentile = models.FloatField(
        verbose_name="SR-HSE_drugbank_approved_percentile", null=True, blank=True
    )
    sr_mmp_drugbank_approved_percentile = models.FloatField(
        verbose_name="SR-MMP_drugbank_approved_percentile", null=True, blank=True
    )
    sr_p53_drugbank_approved_percentile = models.FloatField(
        verbose_name="SR-p53_drugbank_approved_percentile", null=True, blank=True
    )
    skin_reaction_drugbank_approved_percentile = models.FloatField(
        verbose_name="Skin_Reaction_drugbank_approved_percentile", null=True, blank=True
    )
    herg_drugbank_approved_percentile = models.FloatField(
        verbose_name="hERG_drugbank_approved_percentile", null=True, blank=True
    )
    caco2_wang_drugbank_approved_percentile = models.FloatField(
        verbose_name="Caco2_Wang_drugbank_approved_percentile", null=True, blank=True
    )
    clearance_hepatocyte_az_drugbank_approved_percentile = models.FloatField(
        verbose_name="Clearance_Hepatocyte_AZ_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    clearance_microsome_az_drugbank_approved_percentile = models.FloatField(
        verbose_name="Clearance_Microsome_AZ_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    half_life_obach_drugbank_approved_percentile = models.FloatField(
        verbose_name="Half_Life_Obach_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    hydrationfreeenergy_freesolv_drugbank_approved_percentile = models.FloatField(
        verbose_name="HydrationFreeEnergy_FreeSolv_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    ld50_zhu_drugbank_approved_percentile = models.FloatField(
        verbose_name="LD50_Zhu_drugbank_approved_percentile", null=True, blank=True
    )
    lipophilicity_astrazeneca_drugbank_approved_percentile = models.FloatField(
        verbose_name="Lipophilicity_AstraZeneca_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    ppbr_az_drugbank_approved_percentile = models.FloatField(
        verbose_name="PPBR_AZ_drugbank_approved_percentile", null=True, blank=True
    )
    solubility_aqsoldb_drugbank_approved_percentile = models.FloatField(
        verbose_name="Solubility_AqSolDB_drugbank_approved_percentile",
        null=True,
        blank=True,
    )
    vdss_lombardo_drugbank_approved_percentile = models.FloatField(
        verbose_name="VDss_Lombardo_drugbank_approved_percentile", null=True, blank=True
    )

    def save(self, *args, **kwargs):
        if not MoleculeRecord.objects.filter(smiles=self.smiles).exists():
            model = ADMETModel()
            prediction = model.predict(smiles=self.smiles)

            for key, value in prediction.items():
                try:
                    adapted_key = key.lower().replace("-", "_")
                    setattr(self, adapted_key, value)

                except:
                    # We don't want to handle unexpected keys or sth, so just skip errors here
                    pass

        super().save(*args, **kwargs)
