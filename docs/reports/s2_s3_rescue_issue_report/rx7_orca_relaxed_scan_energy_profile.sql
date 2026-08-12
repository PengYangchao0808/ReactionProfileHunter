-- Report-local, deterministic transcription of the optimized-point ORCA
-- relaxscanact ledger for rx7/product_major.  Energies are normalized to the
-- final stretched point in the report snapshot.
WITH energy_profile(point, mean_bond_A, relative_energy_kcal_mol) AS (
  VALUES
    (0, 1.537, -61.454), (1, 1.654, -56.229), (2, 1.771, -44.942),
    (3, 1.888, -32.526), (4, 2.004, -21.989), (5, 2.121, -14.421),
    (6, 2.238, -9.223),  (7, 2.354, -6.070),  (8, 2.471, -4.213),
    (9, 2.588, -3.040),  (10, 2.704, -2.418), (11, 2.821, -2.091),
    (12, 2.938, -1.846), (13, 3.054, -1.705), (14, 3.171, -1.403),
    (15, 3.284, -0.817), (16, 3.400, 0.000)
)
SELECT point, mean_bond_A, relative_energy_kcal_mol
FROM energy_profile
ORDER BY point;
