# coding: utf-8

"""
Lepton selection methods.
"""

from __future__ import annotations

import law

from operator import or_
from functools import reduce

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import (
    set_ak_column, sorted_indices_from_mask, flat_np_view, full_like,
)
from columnflow.util import maybe_import

from multilepton.util import IF_NANO_V9, IF_NANO_GE_V10
from multilepton.config.util import Trigger

np = maybe_import("numpy")
ak = maybe_import("awkward")


logger = law.logger.get_logger(__name__)


def trigger_object_matching(
    vectors1: ak.Array,
    vectors2: ak.Array,
    /,
    *,
    threshold: float = 0.5,
    axis: int = 2,
    event_mask: ak.Array | type(Ellipsis) | None = None,
) -> ak.Array:
    """
    Helper to check per object in *vectors1* if there is at least one object in *vectors2* that
    leads to a delta R metric below *threshold*. The final reduction is applied over *axis* of the
    resulting metric table containing the full combinatorics. If an *event_mask* is given, the
    the matching is performed only for those events, but a full object mask with the same shape as
    that of *vectors1* is returned, which all objects set to *False* where not matching was done.
    """
    # handle event masks
    used_event_mask = event_mask is not None and event_mask is not Ellipsis
    event_mask = Ellipsis if event_mask is None else event_mask

    # delta_r for all combinations
    dr = vectors1[event_mask].metric_table(vectors2[event_mask])

    # check per element in vectors1 if there is at least one matching element in vectors2
    any_match = ak.any(dr < threshold, axis=axis)

    # expand to original shape if an event mask was given
    if used_event_mask:
        full_any_match = full_like(vectors1.pt, False, dtype=bool)
        flat_full_any_match = flat_np_view(full_any_match)
        flat_full_any_match[flat_np_view(full_any_match | event_mask)] = flat_np_view(any_match)
        any_match = full_any_match

    return any_match


def update_channel_ids(
    events: ak.Array,
    previous_channel_ids: ak.Array,
    correct_channel_id: int,
    channel_mask: ak.Array,
) -> ak.Array:
    """
    Check if the events in the is_mask can be inside the given channel
    or have already been sorted in another channel before.
    """
    events_not_in_channel = (previous_channel_ids != 0) & (previous_channel_ids != correct_channel_id)
    channel_id_overwrite = events_not_in_channel & channel_mask
    if ak.any(channel_id_overwrite):
        raise ValueError(
            "The channel_ids of some events are being set to two different values. "
            "The first event of this chunk concerned has index",
            ak.where(channel_id_overwrite)[0],
        )
    return ak.where(channel_mask, correct_channel_id, previous_channel_ids)


@selector(
    uses={
        "Electron.{pt,eta,phi,dxy,dz,pfRelIso03_all,seediEtaOriX,seediPhiOriY}",
        IF_NANO_V9("Electron.mvaFall17V2{Iso_WP80,Iso_WP90}"),
        IF_NANO_GE_V10("Electron.{mvaIso_WP80,mvaIso_WP90}"),
    },
    exposed=False,
)
def electron_selection(
    self: Selector,
    events: ak.Array,
    trigger: Trigger,
    **kwargs,
) -> tuple[ak.Array | None, ak.Array]:
    """
    Electron selection returning two sets of masks for default and veto electrons.
    See https://twiki.cern.ch/twiki/bin/view/CMS/EgammaNanoAOD?rev=4
    """
    is_2016 = self.config_inst.campaign.x.year == 2016
    is_2022_post = (
        self.config_inst.campaign.x.year == 2022 and
        self.config_inst.campaign.has_tag("postEE")
    )
    is_single = trigger.has_tag("single_e")
    is_cross = trigger.has_tag("cross_e_tau")

    # obtain mva flags, which might be located at different routes, depending on the nano version
    if "mvaIso_WP80" in events.Electron.fields:
        # >= nano v10
        # beware that the available Iso should be mvaFall17V2 for run2 files, not Winter22V1,
        # check this in original root files if necessary
        mva_iso_wp80 = events.Electron.mvaIso_WP80
        mva_iso_wp90 = events.Electron.mvaIso_WP90
    else:
        # <= nano v9
        mva_iso_wp80 = events.Electron.mvaFall17V2Iso_WP80
        mva_iso_wp90 = events.Electron.mvaFall17V2Iso_WP90

    # default electron mask
    analysis_mask = None
    control_mask = None
    if is_single or is_cross or True:  # investigate why trigger dependence on providing masks
        min_pt = 26.0 if is_2016 else (31.0 if is_single else 25.0)
        max_eta = 2.5 if is_single else 2.1
        default_mask = (
            (mva_iso_wp80 == 1) &
            (abs(events.Electron.eta) < max_eta) &
            (abs(events.Electron.dxy) < 0.045) &
            (abs(events.Electron.dz) < 0.2)
        )

        # additional cut in 2022 post-EE
        # see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis?rev=162#From_ECAL_and_EGM
        if is_2022_post:
            default_mask = default_mask & ~(
                (events.Electron.eta > 1.556) &
                (events.Electron.seediEtaOriX < 45) &
                (events.Electron.seediPhiOriY > 72)
            )

        # control mask for the electron selection
        control_mask = default_mask & (events.Electron.pt > 10)
        analysis_mask = default_mask & (events.Electron.pt > min_pt)

    # veto electron mask (must be trigger independent!)
    veto_mask = (
        (mva_iso_wp90 == 1) &
        (abs(events.Electron.eta) < 2.5) &
        (abs(events.Electron.dxy) < 0.045) &
        (abs(events.Electron.dz) < 0.2) &
        (events.Electron.pt > 10.0)
    )

    return analysis_mask, control_mask, veto_mask


@electron_selection.init
def electron_selection_init(self) -> None:
    if self.config_inst.campaign.x.run == 3 and self.config_inst.campaign.x.year == 2022:
        self.shifts |= {
            shift_inst.name for shift_inst in self.config_inst.shifts
            if shift_inst.has_tag(("ees", "eer"))
        }


@selector(
    uses={"{Electron,TrigObj}.{pt,eta,phi}"},
    exposed=False,
)
def electron_trigger_matching(
    self: Selector,
    events: ak.Array,
    trigger: Trigger,
    trigger_fired: ak.Array,
    leg_masks: dict[str, ak.Array],
    **kwargs,
) -> tuple[ak.Array]:
    """
    Electron trigger matching.
    """
    is_single = trigger.has_tag("single_e")
    is_cross = trigger.has_tag("cross_e_tau")

    # catch config errors
    assert is_single or is_cross
    assert trigger.n_legs == len(leg_masks) == (1 if is_single else 2)
    assert abs(trigger.legs["e"].pdg_id) == 11

    return trigger_object_matching(
        events.Electron,
        events.TrigObj[leg_masks["e"]],
        event_mask=trigger_fired,
    )


@selector(
    uses={"Muon.{pt,eta,phi,looseId,mediumId,tightId,pfRelIso04_all,dxy,dz}"},
    exposed=False,
)
def muon_selection(
    self: Selector,
    events: ak.Array,
    trigger: Trigger,
    **kwargs,
) -> tuple[ak.Array | None, ak.Array]:
    """
    Muon selection returning two sets of masks for default and veto muons.

    References:

    - Isolation working point: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=59
    - ID und ISO : https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017?rev=15

    relaxed for multilepton, to be replaced with lepMVA later on
    """
    is_2016 = self.config_inst.campaign.x.year == 2016
    is_single = trigger.has_tag("single_mu")
    is_cross = trigger.has_tag("cross_mu_tau")

    # default muon mask
    analysis_mask = None
    control_mask = None
    if is_single or is_cross or True:  # investigate why trigger dependence on providing masks at all
        if is_2016:
            min_pt = 23.0 if is_single else 20.0
        else:
            min_pt = 26.0 if is_single else 22.0
        eta_cut = 2.4 if is_single else 2.1
        default_mask = (
            (events.Muon.mediumId == 1) &
            (abs(events.Muon.eta) < eta_cut) &
            (abs(events.Muon.dxy) < 0.045) &
            (abs(events.Muon.dz) < 0.2) &
            (events.Muon.pfRelIso04_all < 0.4)
        )
        control_mask = default_mask & (events.Muon.pt > 15)  # at the moment can not go below 15 because of muon SF
        analysis_mask = default_mask & (events.Muon.pt > min_pt)

    # veto muon mask (must be trigger independent!)
    veto_mask = (
        ((events.Muon.looseId == 1) | (events.Muon.mediumId == 1)) &
        (abs(events.Muon.eta) < 2.4) &
        (abs(events.Muon.dxy) < 0.045) &
        (abs(events.Muon.dz) < 0.2) &
        (events.Muon.pfRelIso04_all < 0.4) &
        (events.Muon.pt > 10)
    )

    return analysis_mask, control_mask, veto_mask


@selector(
    uses={"{Muon,TrigObj}.{pt,eta,phi}"},
    exposed=False,
)
def muon_trigger_matching(
    self: Selector,
    events: ak.Array,
    trigger: Trigger,
    trigger_fired: ak.Array,
    leg_masks: dict[str, ak.Array],
    **kwargs,
) -> tuple[ak.Array]:
    """
    Muon trigger matching.
    """
    is_single = trigger.has_tag("single_mu")
    is_cross = trigger.has_tag("cross_mu_tau")

    # catch config errors
    assert is_single or is_cross
    assert trigger.n_legs == len(leg_masks) == (1 if is_single else 2)
    assert abs(trigger.legs["mu"].pdg_id) == 13

    return trigger_object_matching(
        events.Muon,
        events.TrigObj[leg_masks["mu"]],
        event_mask=trigger_fired,
    )


@selector(
    uses={
        "Tau.{pt,eta,phi,dz,decayMode}",
        "{Electron,Muon,TrigObj}.{pt,eta,phi}",
    },
    # shifts are declared dynamically below in tau_selection_init
    exposed=False,
)
def tau_selection(
    self: Selector,
    events: ak.Array,
    trigger: Trigger,
    electron_mask: ak.Array | None,
    muon_mask: ak.Array | None,
    **kwargs,
) -> tuple[ak.Array, ak.Array]:
    """
    Tau selection returning a masks for taus that are at least VVLoose isolated (vs jet)
    and a second mask to select isolated ones, eventually to separate normal and iso inverted taus
    for QCD estimations.
    """
    # return empty mask if no tagged taus exists in the chunk
    if ak.all(ak.num(events.Tau) == 0):
        logger.info("no taus found in event chunk")
        false_mask = full_like(events.Tau.pt, False, dtype=bool)
        return false_mask, false_mask

    is_single_e = trigger.has_tag("single_e")
    is_single_mu = trigger.has_tag("single_mu")
    is_cross_e = trigger.has_tag("cross_e_tau")
    is_cross_mu = trigger.has_tag("cross_mu_tau")
    is_cross_tau = trigger.has_tag("cross_tau_tau")
    is_cross_tau_vbf = trigger.has_tag("cross_tau_tau_vbf")
    is_cross_tau_jet = trigger.has_tag("cross_tau_tau_jet")
    is_2016 = self.config_inst.campaign.x.year == 2016
    is_run3 = self.config_inst.campaign.x.run == 3
    get_tau_tagger = lambda tag: f"id{self.config_inst.x.tau_tagger}VS{tag}"
    wp_config = self.config_inst.x.tau_id_working_points

    # determine minimum pt and maximum eta
    max_eta = 2.5
    base_pt = 20.0
    if is_single_e or is_single_mu:
        min_pt = 20.0
    elif is_cross_e:
        # only existing after 2016
        min_pt = 0.0 if is_2016 else 35.0
    elif is_cross_mu:
        min_pt = 25.0 if is_2016 else 32.0
    elif is_cross_tau:
        min_pt = 40.0
    elif is_cross_tau_vbf:
        # only existing after 2016
        min_pt = 0.0 if is_2016 else 25.0
    elif is_cross_tau_jet:
        min_pt = None if not is_run3 else 35.0

    # base tau mask for default and qcd sideband tau
    base_mask = (
        (abs(events.Tau.eta) < max_eta) &
        (events.Tau.pt > base_pt) &
        (abs(events.Tau.dz) < 0.2) &
        reduce(or_, [events.Tau.decayMode == mode for mode in (0, 1, 10, 11)]) &
        (events.Tau[get_tau_tagger("jet")] >= wp_config.tau_vs_jet.vvvloose)
        # vs e and mu cuts are channel dependent and thus applied in the overall lepton selection
    )

    # remove taus with too close spatial separation to previously selected leptons
    if electron_mask is not None:
        base_mask = base_mask & ak.all(events.Tau.metric_table(events.Electron[electron_mask]) > 0.5, axis=2)
    if muon_mask is not None:
        base_mask = base_mask & ak.all(events.Tau.metric_table(events.Muon[muon_mask]) > 0.5, axis=2)

    # trigger dependent cuts
    trigger_specific_mask = base_mask & (events.Tau.pt > min_pt)

    # compute the isolation mask separately as it is used to defined (qcd) categories later on
    iso_mask = events.Tau[get_tau_tagger("jet")] >= wp_config.tau_vs_jet.medium

    return base_mask, trigger_specific_mask, iso_mask


@tau_selection.init
def tau_selection_init(self: Selector) -> None:
    # register tec shifts
    self.shifts |= {
        shift_inst.name
        for shift_inst in self.config_inst.shifts
        if shift_inst.has_tag("tec")
    }

    # Add columns for the right tau tagger
    self.uses |= {
        f"Tau.id{self.config_inst.x.tau_tagger}VS{tag}"
        for tag in ("e", "mu", "jet")
    }


@selector(
    uses={"{Tau,TrigObj}.{pt,eta,phi}"},
    # shifts are declared dynamically below in tau_selection_init
    exposed=False,
)
def tau_trigger_matching(
    self: Selector,
    events: ak.Array,
    trigger: Trigger,
    trigger_fired: ak.Array,
    leg_masks: dict[str, ak.Array],
    **kwargs,
) -> tuple[ak.Array]:
    """
    Tau trigger matching.
    """
    if ak.all(ak.num(events.Tau) == 0):
        logger.info("no taus found in event chunk")
        return full_like(events.Tau.pt, False, dtype=bool)

    is_cross_e = trigger.has_tag("cross_e_tau")
    is_cross_mu = trigger.has_tag("cross_mu_tau")
    is_cross_tau = trigger.has_tag("cross_tau_tau")
    is_cross_tau_vbf = trigger.has_tag("cross_tau_tau_vbf")
    is_cross_tau_jet = trigger.has_tag("cross_tau_tau_jet")
    is_any_cross_tau = is_cross_tau or is_cross_tau_vbf or is_cross_tau_jet
    assert is_cross_e or is_cross_mu or is_any_cross_tau

    # start per-tau mask with trigger object matching per leg
    if is_cross_e or is_cross_mu:
        # catch config errors
        assert trigger.n_legs == len(leg_masks) == 2
        assert abs(trigger.legs["tau"].pdg_id) == 15
        # match leg 1
        return trigger_object_matching(
            events.Tau,
            events.TrigObj[leg_masks["tau"]],
            event_mask=trigger_fired,
        )

    # is_any_cross_tau
    # catch config errors
    assert trigger.n_legs == len(leg_masks) >= 2
    assert abs(trigger.legs["tau1"].pdg_id) == 15
    assert abs(trigger.legs["tau2"].pdg_id) == 15

    # match both legs
    matches_leg0 = trigger_object_matching(
        events.Tau,
        events.TrigObj[leg_masks["tau1"]],
        event_mask=trigger_fired,
    )
    matches_leg1 = trigger_object_matching(
        events.Tau,
        events.TrigObj[leg_masks["tau2"]],
        event_mask=trigger_fired,
    )

    # taus need to be matched to at least one leg, but as a side condition
    # each leg has to have at least one match to a tau
    matches = (
        (matches_leg0 | matches_leg1) &
        ak.any(matches_leg0, axis=1) &
        ak.any(matches_leg1, axis=1)
    )

    return matches


@selector(
    uses={
        electron_selection, electron_trigger_matching, muon_selection, muon_trigger_matching,
        tau_selection, tau_trigger_matching,
        "event", "{Electron,Muon,Tau}.{charge,mass}",
    },
    produces={
        electron_selection, electron_trigger_matching, muon_selection, muon_trigger_matching,
        tau_selection, tau_trigger_matching,
        # new columns
        "channel_id", "leptons_os", "tau2_isolated", "single_triggered", "cross_triggered",
        "matched_trigger_ids",
    },
)
def lepton_selection(
    self: Selector,
    events: ak.Array,
    trigger_results: SelectionResult,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    Combined lepton selection.
    """
    wp_config = self.config_inst.x.tau_id_working_points
    get_tau_tagger = lambda tag: f"id{self.config_inst.x.tau_tagger}VS{tag}"

    # get channels from the config
    print(self.config_inst)
    ch_etau = self.config_inst.get_channel("etau")
    ch_mutau = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")
    ch_ee = self.config_inst.get_channel("ee")
    ch_mumu = self.config_inst.get_channel("mumu")
    ch_emu = self.config_inst.get_channel("emu")
    # new 3l channels
    ch_3e = self.config_inst.get_channel("c3e")
    ch_2emu = self.config_inst.get_channel("c2emu")
    ch_e2mu = self.config_inst.get_channel("ce2mu")
    ch_3mu = self.config_inst.get_channel("c3mu")
    # new 4l channels
    ch_4e = self.config_inst.get_channel("c4e")
    ch_3emu = self.config_inst.get_channel("c3emu")
    ch_2e2mu = self.config_inst.get_channel("c2e2mu")
    ch_e3mu = self.config_inst.get_channel("ce3mu")
    ch_4mu = self.config_inst.get_channel("c4mu")
    # To be implemented
    # #new  3l1tau channels
    # ch_3etau = self.config_inst.get_channel("c3etau")
    # ch_2emutau = self.config_inst.get_channel("c2emutau")
    # ch_e2mutau = self.config_inst.get_channel("ce2mutau")
    # ch_3mutau = self.config_inst.get_channel("c3mutau")
    # #new  2l2tau channels
    # ch_2e2tau = self.config_inst.get_channel("c2e2tau")
    # ch_2mu2tau = self.config_inst.get_channel("c2mu2tau")
    # ch_emu2tau = self.config_inst.get_channel("cemu2tau")
    # # new 1l3tau
    # ch_e3tau = self.config_inst.get_channel("ce3tau")
    # ch_mu3tau = self.config_inst.get_channel("cmu3tau")
    # # new 4tau channel
    # ch_4tau = self.config_inst.get_channel("4tau")

    # prepare vectors for output vectors
    false_mask = (abs(events.event) < 0)
    channel_id = np.uint8(1) * false_mask
    tau2_isolated = false_mask
    leptons_os = false_mask
    single_triggered = false_mask
    cross_triggered = false_mask
    sel_electron_mask = full_like(events.Electron.pt, False, dtype=bool)
    sel_muon_mask = full_like(events.Muon.pt, False, dtype=bool)
    sel_tau_mask = full_like(events.Tau.pt, False, dtype=bool)
    leading_taus = events.Tau[:, :0]
    # store trigger ids for the fully matched triggers
    matched_trigger_ids = []
    # store trigger ids for the matched triggers that do have additional jet legs
    lepton_part_trigger_ids = []

    # indices for sorting taus first by isolation, then by pt
    # for this, combine iso and pt values, e.g. iso 255 and pt 32.3 -> 2550032.3
    f = 10**(np.ceil(np.log10(ak.max(events.Tau.pt))) + 2)
    tau_sorting_key = events.Tau[f"raw{self.config_inst.x.tau_tagger}VSjet"] * f + events.Tau.pt
    tau_sorting_indices = ak.argsort(tau_sorting_key, axis=-1, ascending=False)

    CHANNELS = {
    #id/need/veto/triggers/extrahelpers
    "3e"   : {"id": ch_3e.id,   "need": {"e":3},  "veto": {"mu":0,"tau":0}, "trig": ["single_e"], "helpers": []},
    "4e"   : {"id": ch_4e.id,   "need": {"e":4},  "veto": {"mu":0,"tau":0}, "trig": ["single_e"], "helpers": []},
    "3mu"  : {"id": ch_3mu.id,  "need": {"mu":3}, "veto": {"e":0, "tau":0}, "trig": ["single_mu"],"helpers": []},
    "4mu"  : {"id": ch_4mu.id,  "need": {"mu":4}, "veto": {"e":0, "tau":0}, "trig": ["single_mu"],"helpers": []},

    "2emu" : {"id": ch_2emu.id, "need": {"e":2,"mu":1}, "veto": {}, "trig": ["single_e","single_mu"], "helpers": []},
    "e2mu" : {"id": ch_e2mu.id, "need": {"e":1,"mu":2}, "veto": {}, "trig": ["single_e","single_mu"], "helpers": []},
    "3emu" : {"id": ch_3emu.id, "need": {"e":3,"mu":1}, "veto": {}, "trig": ["single_e","single_mu"], "helpers": []},
    "e3mu" : {"id": ch_e3mu.id, "need": {"e":1,"mu":3}, "veto": {}, "trig": ["single_e","single_mu"], "helpers": []},
    "2e2mu": {"id": ch_2e2mu.id,"need": {"e":2,"mu":2}, "veto": {}, "trig": ["single_e","single_mu"], "helpers": []},
}
 

    # prepare vectors for output vectors
    false_mask = (abs(events.event) < 0)
    channel_id = np.uint8(1) * false_mask
    tau2_isolated = false_mask
    leptons_os = false_mask
    single_triggered = false_mask
    cross_triggered = false_mask
    sel_electron_mask = full_like(events.Electron.pt, False, dtype=bool)
    sel_muon_mask = full_like(events.Muon.pt, False, dtype=bool)
    sel_tau_mask = full_like(events.Tau.pt, False, dtype=bool)
    leading_taus = events.Tau[:, :0]
    
    # indices for sorting taus first by isolation, then by pt
    # for this, combine iso and pt values, e.g. iso 255 and pt 32.3 -> 2550032.3
    f = 10**(np.ceil(np.log10(ak.max(events.Tau.pt))) + 2)
    tau_sorting_key = events.Tau[f"raw{self.config_inst.x.tau_tagger}VSjet"] * f + events.Tau.pt
    tau_sorting_indices = ak.argsort(tau_sorting_key, axis=-1, ascending=False)

    # perform each lepton election step separately per trigger, avoid caching
    #sel_kwargs = {**kwargs, "call_force": True}

#INSERTING THE TWO LOOPS HERE
# ────────────────────────────────────────────────────────────────
# 1 FIRST LOOP – build and cache masks once per fired trigger
# ────────────────────────────────────────────────────────────────

    _trig_cache = {}
    _tid_tags = {}

    e_trig_any  = full_like(events.event, False, dtype=bool)   # we OR all fired flags for single_e here
    mu_trig_any = full_like(events.event, False, dtype=bool)   # we OR all fired flags for single_mu here
    e_match_any = full_like(events.Electron.pt, False, dtype=bool)
    mu_match_any = full_like(events.Muon.pt, False, dtype=bool)
    e_ctrl_single = None; e_mask_single = None; e_veto_single = None
    mu_ctrl_single = None; mu_mask_single = None; mu_veto_single = None


    for trigger, fired, leg_masks in trigger_results.x.trigger_data:

        # Generic selections when is_single==true; Afterwards we need to add the generic selections when is_cross==true

        if trigger.has_tag({"single_e"}) and (e_ctrl_single is None):
            e_mask_single, e_ctrl_single, e_veto_single = self[electron_selection](events, trigger, **kwargs)

        if trigger.has_tag({"single_mu"}) and (mu_ctrl_single is None):
            mu_mask_single, mu_ctrl_single, mu_veto_single = self[muon_selection](events, trigger, **kwargs)

        if not ak.any(fired):
            continue

        e_mask, e_ctrl, e_veto  = self[electron_selection](events, trigger, **kwargs)
        mu_mask, mu_ctrl, mu_veto = self[muon_selection](events, trigger, **kwargs)

        if trigger.has_tag({"single_e"}):
            e_match = self[electron_trigger_matching](events, trigger, fired, leg_masks, **kwargs)
            e_trig_any  = e_trig_any  | fired      # “any single_e fired in this event?”
            e_match_any = e_match_any | e_match    # OR electron matching across all single_e tids
        else:
            # same jagged shape as events.Electron.pt; all False means "no e matched this trigger"
            e_match = full_like(events.Electron.pt, False, dtype=bool)

        # muon matching: only for triggers with a muon leg
        if trigger.has_tag({"single_mu"}):
            mu_match = self[muon_trigger_matching](events, trigger, fired, leg_masks, **kwargs)
            mu_trig_any = mu_trig_any | fired      # “any single_mu fired in this event?”
        else:
            mu_match = full_like(events.Muon.pt, False, dtype=bool)

        tid = trigger.id                                                          #caching information particular to any trigger id
        _trig_cache.update({
            (tid,"e"):e_mask, (tid,"e_ctrl"):e_ctrl, (tid,"e_veto"):e_veto,   
            (tid,"mu"):mu_mask, (tid,"mu_ctrl"):mu_ctrl, (tid,"mu_veto"):mu_veto,
            (tid,"e_match"):e_match, (tid,"mu_match"):mu_match,
        })

        _tid_tags[tid] = set(trigger.tags)

    #Now it is useful to define orthogonal masks: events trigger only on single electrons or single muons 

    e_only        = e_trig_any  & ~mu_trig_any   # only single_e fired
    mu_only       = mu_trig_any & ~e_trig_any    # only single_mu fired
    both_families = e_trig_any  &  mu_trig_any   # both fired

    single_e_tids  = [tid for tid, tags in _tid_tags.items() if "single_e"  in tags]
    single_mu_tids = [tid for tid, tags in _tid_tags.items() if "single_mu" in tags] 

    _trig_cache.update({
        ("fam", "e_trig_any"):  e_trig_any,       #set of events that have triggered at least one single_e trigger
        ("fam", "mu_trig_any"): mu_trig_any,      #set of events that have triggered at least one single_mu trigger
        ("fam", "e_only"):        e_only,         #set of events that have triggered at least one single_e trigger and no one single_mu trigger
        ("fam", "mu_only"):       mu_only,        #set of events that have triggered at least one single_mu trigger and no one single_e trigger
        ("fam", "both_families"): both_families,  #set of events that have triggered at least one singe_e and single_mu trigger

        ("fam", "e_match_any"):   e_match_any,    #Electrons that have matched a single_e trigger object
        
        #the masks below neglect electron/muon selection when is_cross==true (see electron_selection method definition)
        ("fam", "e_ctrl_single"): e_ctrl_single,  #Electrons that passed the control selection provided that is_single==true
        ("fam", "e_mask_single"): e_mask_single,  #Electrons that passed the analysis selection provided that is_single==true
        ("fam", "e_veto_single"): e_veto_single,  #Electrons that passed at veto conditions provided that is_single==true

        ("fam", "mu_ctrl_single"): mu_ctrl_single,
        ("fam", "mu_mask_single"): mu_mask_single,
        ("fam", "mu_veto_single"): mu_veto_single,  
    })
  

# ────────────────────────────────────────────────────────────────
# 2 SECOND LOOP – evaluate every physics channel once
# ────────────────────────────────────────────────────────────────
    for ch_key, spec in CHANNELS.items():

        if ch_key not in {"4e", "4mu", "2e2mu", "3emu", "e3mu"}:
            continue

    # dataset guards identical to bbtautau code
        if ch_key == "4e":
            if not (self.dataset_inst.is_mc or self.dataset_inst.has_tag("ee")):
                continue
            trig_ids = single_e_tids    
        if ch_key == "4mu":
            if not (self.dataset_inst.is_mc or self.dataset_inst.has_tag("mumu")):
                continue
            trig_ids = single_mu_tids
        if ch_key == "2e2mu":
            if self.dataset_inst.has_tag("emu_from_e"):
                trig_ids = single_e_tids
            elif self.dataset_inst.has_tag("emu_from_mu"):
                trig_ids = single_mu_tids
            elif self.dataset_inst.is_mc:
                trig_ids = single_e_tids + single_mu_tids
            else :
                continue
        if ch_key == "3emu":
            if self.dataset_inst.has_tag("emu_from_e"):
                trig_ids = single_e_tids
            elif self.dataset_inst.has_tag("emu_from_mu"):
                trig_ids = single_mu_tids
            elif self.dataset_inst.is_mc:
                trig_ids = single_e_tids + single_mu_tids
            else:
                continue
        if ch_key == "e3mu":
            if self.dataset_inst.has_tag("emu_from_e"):
                trig_ids = single_e_tids
            elif self.dataset_inst.has_tag("emu_from_mu"):
                trig_ids = single_mu_tids
            elif self.dataset_inst.is_mc:
                trig_ids = single_e_tids + single_mu_tids
            else:
                continue

        good_evt = ak.zeros_like(events.event, dtype=bool)

        for tid in trig_ids:
            e_mask  = _trig_cache[(tid,"e")];   e_ctrl  = _trig_cache[(tid,"e_ctrl")]
            mu_mask = _trig_cache[(tid,"mu")];  mu_ctrl = _trig_cache[(tid,"mu_ctrl")]
            e_veto  = _trig_cache[(tid,"e_veto")]; mu_veto = _trig_cache[(tid,"mu_veto")]
            e_match = _trig_cache[(tid,"e_match")]; mu_match = _trig_cache[(tid,"mu_match")]

            ok = ak.ones_like(events.event, dtype=bool)
              
            if ch_key == "4e":
                base_ok = (
                    (ak.sum(e_mask,  axis=1) >= 1) &
                    (ak.sum(e_ctrl,  axis=1) == 4) &
                    (ak.sum(e_veto,  axis=1) == 4) &   
                    (ak.sum(mu_veto, axis=1) == 0) &
                    ak.any(e_match & e_mask, axis=1)
                )

                ok = ak.where(base_ok, ok, False)

                for flav, maxn in spec["veto"].items():
                    if flav not in {"e", "mu"}:
                        continue
                    veto_mask = {"e": e_veto, "mu": mu_veto}[flav]
                    ok = ak.where(ak.sum(veto_mask, axis=1) <= maxn, ok, False)

                charge = events.Electron[e_ctrl & ok[:, None]].charge
                os = (ak.num(charge) == 4) & (ak.sum(charge, axis=1) == 0)
                os = ak.fill_none(os, False)

                leptons_os       = ak.where(ok, leptons_os | os, leptons_os)
                single_triggered = ak.where(ok, True, single_triggered)
                sel_electron_mask = ak.where(ok, sel_electron_mask | e_ctrl | e_mask, sel_electron_mask)


            elif ch_key == "4mu":
                base_ok = (
                    (ak.sum(mu_mask,  axis=1) >= 1) &               # at least one analysis muon present
                    (ak.sum(mu_ctrl,  axis=1) == 4) &               # exactly four control muons
                    (ak.sum(mu_veto,  axis=1) == 4) &               # exactly four veto muons
                    (ak.sum(e_veto,   axis=1) == 0) &               # zero veto electrons
                    ak.any(mu_match & mu_mask, axis=1)              # trigger matching with the analysis muon
                )

                ok = ak.where(base_ok, ok, False)

    # additional category vetoes
                for flav, maxn in spec["veto"].items():
                    if flav not in {"e", "mu"}:
                        continue
                    veto_mask = {"e": e_veto, "mu": mu_veto}[flav]
                    ok = ak.where(ak.sum(veto_mask, axis=1) <= maxn, ok, False)

    # opposite-sign
                charge = events.Muon[mu_ctrl & ok[:, None]].charge
                os = (ak.num(charge) == 4) & (ak.sum(charge, axis=1) == 0)
                os = ak.fill_none(os, False)

    # final objects
                leptons_os      = ak.where(ok, leptons_os | os, leptons_os)
                single_triggered = ak.where(ok, True, single_triggered)
                sel_muon_mask    = ak.where(ok, sel_muon_mask | mu_ctrl | mu_mask, sel_muon_mask)

            elif ch_key == "2e2mu":

                if tid in single_e_tids:
    # emu_from_e — accept ONLY events with e_only (anti-overlap)
                    ok = ak.where(e_only, ok, False)

                    trig_electron_mask = e_mask & e_match

    # exactly 2 control per flavor, exactly 2 veto per flavor
                    ok = ak.where(ak.sum(e_ctrl,           axis=1) == 2, ok, False)
                    ok = ak.where(ak.sum(e_veto,           axis=1) == 2, ok, False)
                    ok = ak.where(ak.sum(mu_ctrl,   axis=1) == 2, ok, False)
                    ok = ak.where(ak.sum(mu_veto,   axis=1) == 2, ok, False)

    # trigger-side: >=1 analysis e and >=1 matched e (this tid)
                    ok = ak.where(ak.sum(e_mask,           axis=1) >= 1, ok, False)
                    ok = ak.where(ak.sum(trig_electron_mask, axis=1) >= 1, ok, False)

    # muons: at least one offline 
                    ok = ak.where(ak.sum(mu_ctrl_single,   axis=1) >= 1, ok, False)

                    single_triggered  = ak.where(ok, True, single_triggered)
                    sel_electron_mask = ak.where(ok, sel_electron_mask | e_ctrl,         sel_electron_mask)
                    sel_muon_mask     = ak.where(ok, sel_muon_mask     | mu_ctrl_single, sel_muon_mask)

                elif tid in single_mu_tids:
    # emu_from_mu: allow both_families, the matching/logic below handles e-side as well
                    trig_muon_mask = mu_mask & mu_match

                    ok = ak.where(ak.sum(mu_ctrl, axis=1) == 2, ok, False)
                    ok = ak.where(ak.sum(mu_veto, axis=1) == 2, ok, False)
                    ok = ak.where(ak.sum(e_ctrl,  axis=1) == 2, ok, False)
                    ok = ak.where(ak.sum(e_veto,  axis=1) == 2, ok, False)

                    ok = ak.where(ak.sum(mu_mask,         axis=1) >= 1, ok, False)
                    ok = ak.where(ak.sum(trig_muon_mask,  axis=1) >= 1, ok, False)

    # electrons: if any single-e fired, require analysis+matching, else control is enough
                    emu_electron_mask  = ak.where(e_trig_any, e_mask_single, e_ctrl_single)
                    e_match_mask       = (e_match_any | ~e_trig_any)
                    trig_electron_mask = emu_electron_mask & e_match_mask

                    ok = ak.where(ak.sum(emu_electron_mask,  axis=1) >= 1, ok, False)
                    ok = ak.where(ak.sum(trig_electron_mask, axis=1) >= 1, ok, False)

                    single_triggered  = ak.where(ok, True, single_triggered)
                    sel_electron_mask = ak.where(ok, sel_electron_mask | e_ctrl_single, sel_electron_mask)
                    sel_muon_mask     = ak.where(ok, sel_muon_mask     | mu_ctrl,       sel_muon_mask)
                

            elif ch_key == "3emu":

                if tid in single_e_tids:
    # emu_from_e — accept ONLY events with e_only (anti-overlap)
                    ok = ak.where(e_only, ok, False)

                    trig_electron_mask = e_mask & e_match

                    ok = ak.where(ak.sum(e_ctrl,           axis=1) == 3, ok, False)
                    ok = ak.where(ak.sum(e_veto,           axis=1) == 3, ok, False)
                    ok = ak.where(ak.sum(mu_ctrl,   axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(mu_veto,   axis=1) == 1, ok, False)

    # trigger-side: >=1 analysis e and >=1 matched e (this tid)
                    ok = ak.where(ak.sum(e_mask,           axis=1) >= 1, ok, False)
                    ok = ak.where(ak.sum(trig_electron_mask, axis=1) >= 1, ok, False)

    # muons: at least one offline
                    ok = ak.where(ak.sum(mu_ctrl_single,   axis=1) == 1, ok, False)

                    single_triggered  = ak.where(ok, True, single_triggered)
                    sel_electron_mask = ak.where(ok, sel_electron_mask | e_ctrl,         sel_electron_mask)
                    sel_muon_mask     = ak.where(ok, sel_muon_mask     | mu_ctrl_single, sel_muon_mask)

                elif tid in single_mu_tids:
    # emu_from_mu — allow both_families; the matching/logic below handles e-side
                    trig_muon_mask = mu_mask & mu_match

                    ok = ak.where(ak.sum(mu_ctrl, axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(mu_veto, axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(e_ctrl,  axis=1) == 3, ok, False)
                    ok = ak.where(ak.sum(e_veto,  axis=1) == 3, ok, False)

                    ok = ak.where(ak.sum(mu_mask,         axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(trig_muon_mask,  axis=1) == 1, ok, False)

    # electrons: if any single-e fired, require analysis+matching; else control is enough
                    emu_electron_mask  = ak.where(e_trig_any, e_mask_single, e_ctrl_single)
                    e_match_mask       = (e_match_any | ~e_trig_any)
                    trig_electron_mask = emu_electron_mask & e_match_mask

                    ok = ak.where(ak.sum(emu_electron_mask,  axis=1) >= 1, ok, False)
                    ok = ak.where(ak.sum(trig_electron_mask, axis=1) >= 1, ok, False)

                    single_triggered  = ak.where(ok, True, single_triggered)
                    sel_electron_mask = ak.where(ok, sel_electron_mask | e_ctrl_single, sel_electron_mask)
                    sel_muon_mask     = ak.where(ok, sel_muon_mask     | mu_ctrl,       sel_muon_mask)


            elif ch_key == "e3mu":

                if tid in single_e_tids:
    # emu_from_e — accept ONLY events with e_only (anti-overlap)
                    ok = ak.where(e_only, ok, False)

                    trig_electron_mask = e_mask & e_match

                    ok = ak.where(ak.sum(e_ctrl,           axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(e_veto,           axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(mu_ctrl,   axis=1) == 3, ok, False)
                    ok = ak.where(ak.sum(mu_veto,   axis=1) == 3, ok, False)

                    ok = ak.where(ak.sum(e_mask,           axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(trig_electron_mask, axis=1) == 1, ok, False)

    # muons: at least one offline (
                    ok = ak.where(ak.sum(mu_ctrl_single,   axis=1) >= 1, ok, False)

                    single_triggered  = ak.where(ok, True, single_triggered)
                    sel_electron_mask = ak.where(ok, sel_electron_mask | e_ctrl,         sel_electron_mask)
                    sel_muon_mask     = ak.where(ok, sel_muon_mask     | mu_ctrl_single, sel_muon_mask)

                elif tid in single_mu_tids:
    # emu_from_mu — allow both_families; the matching/logic below handles e-side
                    trig_muon_mask = mu_mask & mu_match

                    ok = ak.where(ak.sum(mu_ctrl, axis=1) == 3, ok, False)
                    ok = ak.where(ak.sum(mu_veto, axis=1) == 3, ok, False)
                    ok = ak.where(ak.sum(e_ctrl,  axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(e_veto,  axis=1) == 1, ok, False)

                    ok = ak.where(ak.sum(mu_mask,         axis=1) >= 1, ok, False)
                    ok = ak.where(ak.sum(trig_muon_mask,  axis=1) >= 1, ok, False)

    # electrons: if any single-e fired, require analysis+matching; else control is enough
                    emu_electron_mask  = ak.where(e_trig_any, e_mask_single, e_ctrl_single)
                    e_match_mask       = (e_match_any | ~e_trig_any)
                    trig_electron_mask = emu_electron_mask & e_match_mask

                    ok = ak.where(ak.sum(emu_electron_mask,  axis=1) == 1, ok, False)
                    ok = ak.where(ak.sum(trig_electron_mask, axis=1) == 1, ok, False)

                    single_triggered  = ak.where(ok, True, single_triggered)
                    sel_electron_mask = ak.where(ok, sel_electron_mask | e_ctrl_single, sel_electron_mask)
                    sel_muon_mask     = ak.where(ok, sel_muon_mask     | mu_ctrl,       sel_muon_mask)                    
                

        # accumulate over triggers
            good_evt = ak.where(ok, True, good_evt)

        channel_id = update_channel_ids(events, channel_id, spec["id"], good_evt)



    # some final type conversions
    channel_id = ak.values_astype(channel_id, np.uint8)
    leptons_os = ak.fill_none(leptons_os, False)

    # concatenate matched trigger ids
    empty_ids = ak.singletons(full_like(events.event, 0, dtype=np.int32), axis=0)[:, :0]
    merge_ids = lambda ids: ak.values_astype(ak.concatenate(ids, axis=1), np.int32) if ids else empty_ids
    matched_trigger_ids = merge_ids(matched_trigger_ids)
    lepton_part_trigger_ids = merge_ids(lepton_part_trigger_ids)

    # save new columns
    events = set_ak_column(events, "channel_id", channel_id)
    events = set_ak_column(events, "leptons_os", leptons_os)
    events = set_ak_column(events, "tau2_isolated", tau2_isolated)
    events = set_ak_column(events, "single_triggered", single_triggered)
    events = set_ak_column(events, "cross_triggered", cross_triggered)
    events = set_ak_column(events, "matched_trigger_ids", matched_trigger_ids)

    # convert lepton masks to sorted indices (pt for e/mu, iso for tau)
    sel_electron_indices = sorted_indices_from_mask(sel_electron_mask, events.Electron.pt, ascending=False)
    sel_muon_indices = sorted_indices_from_mask(sel_muon_mask, events.Muon.pt, ascending=False)
    sel_tau_indices = sorted_indices_from_mask(sel_tau_mask, tau_sorting_key, ascending=False)

    return events, SelectionResult(
        steps={
            "lepton": channel_id != 0,
        },
        objects={
            "Electron": {
                "Electron": sel_electron_indices,
            },
            "Muon": {
                "Muon": sel_muon_indices,
            },
            "Tau": {
                "Tau": sel_tau_indices,
            },
        },
        aux={
            # save the selected lepton pair for the duration of the selection
            # multiplication of a coffea particle with 1 yields the lorentz vector
            "lepton_pair": ak.concatenate(
                [
                    events.Electron[sel_electron_indices] * 1,
                    events.Muon[sel_muon_indices] * 1,
                    events.Tau[sel_tau_indices] * 1,
                ],
                axis=1,
            )[:, :2],

            # save the matched trigger ids of the trigger with jet legs for the duration of the selection
            # these will be updated in the jet selection and then stored in the matched_trigger_ids column
            "lepton_part_trigger_ids": lepton_part_trigger_ids,


            # save the leading taus for the duration of the selection
            # exactly 1 for etau/mutau and exactly 2 for tautau
            "leading_taus": leading_taus,
        },
    )


@lepton_selection.init
def lepton_selection_init(self: Selector, **kwargs) -> None:
    # add column to load the raw tau tagger score
    self.uses.add(f"Tau.raw{self.config_inst.x.tau_tagger}VSjet")
