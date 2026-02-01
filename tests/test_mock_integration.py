"""Mock integration tests (lightweight, no external QC).

These cover basic data plumbing and simple helpers without invoking Gaussian/ORCA.
"""


def test_spmatrixreport_activation_energy_prefers_gibbs() -> None:
    from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport

    report = SPMatrixReport(
        g_ts=10.0,
        g_reactant=-10.0,
        g_product=-5.0,
        e_ts_final=0.0,
        e_reactant=0.0,
        e_product=0.0,
    )

    assert report.get_activation_energy() == 20.0


def test_spmatrixreport_reaction_energy_prefers_gibbs() -> None:
    from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport

    report = SPMatrixReport(
        g_ts=0.0,
        g_reactant=-10.0,
        g_product=-5.0,
        e_ts_final=0.0,
        e_reactant=0.0,
        e_product=0.0,
    )

    assert report.get_reaction_energy() == 5.0
