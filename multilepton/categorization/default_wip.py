# coding: utf-8

"""
HH -> multi-leptons selection methods.
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import
import inspect

ak = maybe_import("awkward")

MultiLeptonsChannels = [
    "etau", "mutau", "tautau", "ee", "mumu", "emu",
    "3e", "2emu", "e2mu", "3mu", "4e", "3emu", "2e2mu", "e3mu", "4mu",
    "3etau", "2emutau", "e2mutau", "3mutau",
    "2e2tau", "emu2tau", "2mu2tau", "e3tau", "mu3tau", "4tau",
    "2ess", "emuss", "2muss"
]
# exceptions to the "c{name}" rule (if any don't follow the convention)
CHANNEL_EXCEPTIONS = {
    # Example: "mutau" : "cmutau_alt"
}


@categorizer(uses={"event"})
def cat_all(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    """Keep all events."""
    return events, ak.ones_like(events.event, dtype=bool)

def _make_channel_categorizer(name: str, channel_key: str):
    @categorizer(uses={"channel_id"})
    def func(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
        if not hasattr(self.config_inst, "channels") or not hasattr(self.config_inst.channels.n, channel_key):
            # fallback: mark all as False until config is ready
            return events, ak.zeros_like(events.channel_id, dtype=bool)
        return events, events.channel_id == getattr(self.config_inst.channels.n, channel_key).id

    func.__name__ = f"cat_{name}"
    func.__qualname__ = f"cat_{name}"
    func.__doc__ = f"Select events belonging to the {name} channel."
    return func

def _register_channel_categorizers():
    """Scan config_inst.channels.n and generate all matching categorizer functions."""
    # Loop over attributes in channels.n
    for attr in dir(Categorizer.config_inst.channels.n):  # type: ignore
        if attr.startswith("c"):  # only consider channel-like entries
            name = attr[1:]  # e.g. cetau -> etau
            # Skip special/internal attributes
            if not hasattr(Categorizer.config_inst.channels.n, attr):
                continue
            func = _make_channel_categorizer(name, attr)
            globals()[func.__name__] = func

def register_multilepton_categorizers():
    try:
        _register_channel_categorizers()
    except Exception as e:
        print(f"[INFO] Could not register categorizers yet: {e}")

# ------------------------------------------------------------------------
# FALLBACK: If config_inst is not yet initialized (e.g. during import)
# just define the functions based on known channel names
# ------------------------------------------------------------------------
# Generate these upfront so the file works standalone
for name in MultiLeptonsChannels:
    chkey = CHANNEL_EXCEPTIONS.get(name, f"c{name}")
    globals()[f"cat_{name}"] = _make_channel_categorizer(name, chkey)

register_multilepton_categorizers()
# ------------------------------------------------------------------------
# other Categories 
# ------------------------------------------------------------------------
# 3l/4l inclusive, later split into CR / SR via Z-peak
@categorizer(uses={"channel_id"})
def cat_3l0tau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    catmask = events.channel_id == self.config_inst.channels.n.c3e.id
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c3mu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c2emu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.ce2mu.id)
    return events, catmask

@categorizer(uses={"channel_id"})
def cat_4l(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    catmask = events.channel_id == self.config_inst.channels.n.c4e.id
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c3emu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c2e2mu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.ce3mu.id)
    catmask = catmask | (events.channel_id == self.config_inst.channels.n.c4mu.id)
    return events, catmask

# bveto
@categorizer(uses={"Jet.btagPNetB"})
def cat_bveto_on(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    wp_loose = self.config_inst.x.btag_working_points["particleNet"]["loose"]
    wp_medium = self.config_inst.x.btag_working_points["particleNet"]["medium"]
    tagged_loose = events.Jet.btagPNetB > wp_loose
    tagged_medium = events.Jet.btagPNetB > wp_medium
    veto = (ak.sum(tagged_loose, axis=1) < 2) & (ak.sum(tagged_medium, axis=1) < 1)
    return events, veto

@categorizer(uses={"Jet.btagPNetB"})
def cat_bveto_off(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    wp_loose = self.config_inst.x.btag_working_points["particleNet"]["loose"]
    wp_medium = self.config_inst.x.btag_working_points["particleNet"]["medium"]
    tagged_loose = events.Jet.btagPNetB > wp_loose
    tagged_medium = events.Jet.btagPNetB > wp_medium
    nonveto = (ak.sum(tagged_loose, axis=1) >= 2) | (ak.sum(tagged_medium, axis=1) >= 1)
    return events, nonveto

# The BDT category overlaps with our channels, so we need tight/trigger-matched flags individual for this cat
@categorizer(uses={"ok_bdt_eormu"})
def cat_e_or_mu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.ok_bdt_eormu == 1

@categorizer(uses={"tight_sel_bdt"})
def cat_tight_bdt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # tight true
    return events, events.tight_sel_bdt == 1

@categorizer(uses={"tight_sel_bdt"})
def cat_nontight_bdt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # tight false
    return events, events.tight_sel_bdt == 0

@categorizer(uses={"trig_match_bdt"})
def cat_trigmatch_bdt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # trig match
    return events, events.trig_match_bdt == 1

@categorizer(uses={"trig_match_bdt"})
def cat_nontrigmatch_bdt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # trig match false
    return events, events.trig_match_bdt == 0

# Tight and trigger matching flags for the physical channels
@categorizer(uses={"tight_sel"})
def cat_tight(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # tight true
    return events, events.tight_sel == 1

@categorizer(uses={"tight_sel"})
def cat_nontight(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # tight false
    return events, events.tight_sel == 0

@categorizer(uses={"trig_match"})
def cat_trigmatch(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # trig match
    return events, events.trig_match == 1

@categorizer(uses={"trig_match"})
def cat_nontrigmatch(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # trig match false
    return events, events.trig_match == 0

# QCD regions
@categorizer(uses={"leptons_os"})
def cat_os(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # oppositive sign leptons
    return events, events.leptons_os == 1

@categorizer(uses={"leptons_os"})
def cat_ss(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # same sign leptons
    return events, events.leptons_os == 0

@categorizer(uses={"tau2_isolated"})
def cat_iso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # isolated tau2
    return events, events.tau2_isolated == 1

@categorizer(uses={"tau2_isolated"})
def cat_noniso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # noon-isolated tau2
    return events, events.tau2_isolated == 0

# kinematic regions
@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1

@categorizer(uses={"Jet.{pt,phi}"})
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # two or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 2

@categorizer(uses={"Jet.btagPNetB"})
def cat_res1b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # exactly pnet b-tags
    wp = self.config_inst.x.btag_working_points["particleNet"]["medium"]
    tagged = events.Jet.btagPNetB > wp
    return events, ak.sum(tagged, axis=1) == 1

@categorizer(uses={"Jet.btagPNetB"})
def cat_res2b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # at least two medium pnet b-tags
    wp = self.config_inst.x.btag_working_points["particleNet"]["medium"]
    tagged = events.Jet.btagPNetB > wp
    return events, ak.sum(tagged, axis=1) >= 2

@categorizer(uses={cat_res1b, cat_res2b, "FatJet.{pt,phi}"})
def cat_boosted(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # not res1b or res2b, and exactly one selected fat jet that should also pass a tighter pt cut
    # note: this is just a draft
    mask = (
        (ak.num(events.FatJet, axis=1) == 1) &
        (ak.sum(events.FatJet.pt > 350, axis=1) == 1) &
        ~self[cat_res1b](events, **kwargs)[1] &
        ~self[cat_res2b](events, **kwargs)[1]
    )
    return events, mask

@categorizer(uses={"{Electron,Muon,Tau}.{pt,eta,phi,mass}"})
def cat_dy(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # e/mu driven DY region: mll > 40 and met < 30 (to supress tau decays into e/mu)
    leps = ak.concatenate([events.Electron * 1, events.Muon * 1, events.Tau * 1], axis=1)[:, :2]
    mask = (
        (leps.sum(axis=1).mass > 40) &
        (events[self.config_inst.x.met_name].pt < 30)
    )
    return events, mask

@cat_dy.init
def cat_dy_init(self: Categorizer) -> None:
    self.uses.add(f"{self.config_inst.x.met_name}.{{pt,phi}}")

@categorizer(uses={"{Electron,Muon,Tau}.{pt,eta,phi,mass}"})
def cat_tt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # tt region: met > 30 (due to neutrino presence in leptonic w decays)
    mask = events[self.config_inst.x.met_name].pt > 30
    return events, mask

@cat_tt.init
def cat_tt_init(self: Categorizer) -> None:
    self.uses.add(f"{self.config_inst.x.met_name}.{{pt,phi}}")


