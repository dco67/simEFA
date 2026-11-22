"""Captures d'écran de simEFA pour le manuel.
Usage : python3 captures/captures.py --url http://127.0.0.1:8841 --images images [--fichier cap-x.png]
L'application doit tourner (shiny::runApp). Toutes les captures partent de la graine 1234."""
import argparse, sys, time
from pathlib import Path
from playwright.sync_api import sync_playwright

def attendre(page, pause=0.8):
    page.wait_for_function("() => window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected() && !document.documentElement.classList.contains('shiny-busy')", timeout=300000)
    time.sleep(pause)
    page.wait_for_function("() => !document.documentElement.classList.contains('shiny-busy')", timeout=300000)

def onglet(page, titre):
    loc = page.locator(f"xpath=//a[(@data-bs-toggle='tab' or @data-toggle='tab') and normalize-space()=\"{titre}\"]")
    for k in range(loc.count()):
        if loc.nth(k).is_visible():
            loc.nth(k).click(); attendre(page, 1.5); return
    raise RuntimeError("onglet introuvable : " + titre)

def js(page, code, arg=None):
    page.evaluate(code, arg)

def cliquer(page, i, pause=1.5):
    js(page, "id => document.getElementById(id).click()", i); attendre(page, pause)

def saisir(page, i, v):
    js(page, "([id, v]) => { const el = document.getElementById(id); el.value = v; $(el).trigger('change'); }", [i, str(v)]); attendre(page, 0.6)

def choisir(page, i, v):
    js(page, "([id, v]) => { const el=$('#'+id); if (el[0].selectize) el[0].selectize.setValue(v); else el.val(v).trigger('change'); }", [i, v]); attendre(page, 1.0)

def radio(page, nom, v):
    js(page, "([n, v]) => { const el = document.querySelector(`input[name='${n}'][value='${v}']`); el.checked = true; $(el).trigger('change'); }", [nom, v]); attendre(page, 0.8)

def nettoyer(page):
    js(page, "() => document.querySelectorAll('.shiny-notification').forEach(n => n.remove())")

def carte_de(page, ul_id):
    """La carte (card) qui contient la barre d'onglets ul#ul_id."""
    return page.locator(f"#{ul_id}").locator("xpath=ancestor::div[contains(concat(' ', normalize-space(@class), ' '), ' card ')][1]")

def capturer(loc, fichier):
    loc.scroll_into_view_if_needed(); time.sleep(0.4)
    loc.screenshot(path=str(fichier)); print("  ✓", Path(fichier).name, flush=True)

# --- Étapes communes -----------------------------------------------------
def generer(page, scenario="clear", **reglages):
    if scenario != "clear":
        choisir(page, "sim_scenario", scenario)
    for k, v in reglages.items():
        saisir(page, k, v)
    cliquer(page, "generate_example", 2.5); nettoyer(page)

def classique(page, nom):
    radio(page, "import_method", "classic"); choisir(page, "classic_dataset", nom)
    cliquer(page, "load_classic", 2.5); nettoyer(page)

def afe(page, q=None, rotation=None):
    if q is not None: saisir(page, "efa_nfactors", q)
    if rotation is not None: choisir(page, "dynamic_rotation", rotation)
    cliquer(page, "run_efa", 4); nettoyer(page)

# --- Les captures --------------------------------------------------------
def cap_interface(page, out):
    page.screenshot(path=str(out / "cap-interface.png")); print("  ✓ cap-interface.png")

def cap_correlation(page, out):
    generer(page); onglet(page, "Aperçu des données"); onglet(page, "Statistiques descriptives")
    onglet(page, "Matrice de corrélation"); capturer(carte_de(page, "stats_tabs"), out / "cap-correlation.png")

def cap_preparation(page, out):
    classique(page, "rse")
    js(page, "() => { const b = document.querySelector(\".accordion-item[data-value='item_preparation'] .accordion-button\"); if (b.classList.contains('collapsed')) b.click(); }"); time.sleep(0.8)
    choisir(page, "reverse_items", ["Q3", "Q5", "Q8", "Q9", "Q10"])
    cliquer(page, "reverse_selected_items", 2); nettoyer(page)
    capturer(page.locator(".accordion-item[data-value='item_preparation']"), out / "cap-preparation.png")

def cap_prelim(page, out):
    generer(page); onglet(page, "Analyse factorielle"); onglet(page, "Tests préliminaires"); attendre(page, 3)
    capturer(carte_de(page, "efa_tabs"), out / "cap-prelim.png")

def cap_comparaison(page, out):
    generer(page); afe(page); onglet(page, "Analyse factorielle"); onglet(page, "Comparaison des modèles")
    capturer(carte_de(page, "efa_tabs"), out / "cap-comparaison.png")

def cap_charges(page, out):
    generer(page); afe(page); onglet(page, "Interprétation"); onglet(page, "1. Charges")
    capturer(carte_de(page, "interpretation_tabs"), out / "cap-charges.png")

def cap_qualite(page, out):
    generer(page); afe(page); onglet(page, "Interprétation"); onglet(page, "4. Qualité de la solution"); attendre(page, 2)
    capturer(carte_de(page, "interpretation_tabs"), out / "cap-qualite.png")

def cap_structure(page, out):
    classique(page, "hs1939"); afe(page, 3); onglet(page, "Interprétation"); onglet(page, "2. Structure des items")
    capturer(carte_de(page, "interpretation_tabs"), out / "cap-structure.png")

def cap_correles(page, out):
    generer(page, "correlated"); afe(page); onglet(page, "Interprétation"); onglet(page, "3. Facteurs corrélés")
    capturer(carte_de(page, "interpretation_tabs"), out / "cap-correles.png")

def cap_fiabilite(page, out):
    generer(page); afe(page); onglet(page, "Analyse factorielle"); onglet(page, "Fiabilité des facteurs"); attendre(page, 2)
    capturer(carte_de(page, "efa_tabs"), out / "cap-fiabilite.png")

def cap_bootstrap(page, out):
    generer(page); afe(page); onglet(page, "Bootstrap"); saisir(page, "bootstrap_samples", 200)
    js(page, "() => { Math.seedrandom && Math.seedrandom(1); }")
    cliquer(page, "run_bootstrap", 6); nettoyer(page)
    choisir(page, "bootstrap_stability_var", "X1"); attendre(page, 2)
    capturer(page.locator("#bootstrap_loadings_table").locator("xpath=ancestor::div[contains(concat(' ', normalize-space(@class), ' '), ' bslib-grid ')][1]"), out / "cap-bootstrap.png")

def cap_cfa(page, out):
    generer(page); afe(page); onglet(page, "Validation CFA")
    cliquer(page, "generate_validation_sample", 2); cliquer(page, "generate_cfa", 1.5); nettoyer(page)
    onglet(page, "Modèle lavaan")
    capturer(carte_de(page, "cfa_main_tabs"), out / "cap-cfa-syntaxe.png")
    cliquer(page, "run_cfa", 4); nettoyer(page)
    onglet(page, "Résultats et interprétation"); onglet(page, "Charges standardisées"); attendre(page, 2)
    capturer(carte_de(page, "cfa_main_tabs"), out / "cap-cfa-resultats.png")

CAPTURES = {
    "cap-interface.png": cap_interface, "cap-correlation.png": cap_correlation,
    "cap-preparation.png": cap_preparation, "cap-prelim.png": cap_prelim,
    "cap-comparaison.png": cap_comparaison, "cap-charges.png": cap_charges,
    "cap-qualite.png": cap_qualite, "cap-structure.png": cap_structure,
    "cap-correles.png": cap_correles, "cap-fiabilite.png": cap_fiabilite,
    "cap-bootstrap.png": cap_bootstrap, "cap-cfa-syntaxe.png": cap_cfa,
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True); ap.add_argument("--images", default="images")
    ap.add_argument("--fichier"); ap.add_argument("--chrome", default="/opt/pw-browsers/chromium-1194/chrome-linux/chrome")
    a = ap.parse_args(); out = Path(a.images); out.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as p:
        nav = p.chromium.launch(executable_path=a.chrome) if Path(a.chrome).exists() else p.chromium.launch()
        for nom, f in CAPTURES.items():
            if a.fichier and nom != a.fichier: continue
            hauteur = 900 if nom == "cap-interface.png" else 2600
            page = nav.new_page(viewport={"width": 1440, "height": hauteur}, device_scale_factor=1.5, locale="fr-CA")
            page.goto(a.url); attendre(page, 2.5)
            f(page, out); page.close()
        nav.close()

if __name__ == "__main__":
    main()
