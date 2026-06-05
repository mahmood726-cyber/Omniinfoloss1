// tests.js — pure-Node test harness for engine.js (786-M21).
// Expected values are hand-derived independently (see comments), NOT read back
// from the implementation. Run:  node tests.js
const E = require('./engine.js');
const { Stat, Matrix, InfoMetrics, DistributionRecovery, Pooling } = E;

let passed = 0, failed = 0;
function ok(name, cond) { if (cond) { passed++; } else { failed++; console.log('FAIL: ' + name); } }
function approx(name, got, exp, tol = 1e-6) {
    const good = Number.isFinite(got) && Math.abs(got - exp) <= tol;
    if (good) passed++; else { failed++; console.log(`FAIL: ${name}  got=${got} expected=${exp} tol=${tol}`); }
}

// ---------------------------------------------------------------------------
// 1. normalCDF reference points. Stat.pnorm(z)=0.5*(1+erf(z/sqrt2)).
//    Φ(0) MUST be 0.5 (a 0 here is the classic erf/CDF-structure bug). Φ(1.96)≈0.975.
// ---------------------------------------------------------------------------
approx('pnorm(0)=0.5', Stat.pnorm(0), 0.5, 1e-9);
approx('pnorm(1.96)~0.975', Stat.pnorm(1.96), 0.975, 1e-3);
approx('pnorm(-1.96)~0.025', Stat.pnorm(-1.96), 0.025, 1e-3);
ok('pnorm monotone', Stat.pnorm(1) > Stat.pnorm(0) && Stat.pnorm(0) > Stat.pnorm(-1));
// qnorm is the inverse of pnorm
approx('qnorm(0.975)~1.96', Stat.qnorm(0.975), 1.959964, 1e-3);
approx('invLogit(0)=0.5', Stat.invLogit(0), 0.5, 1e-12);

// ---------------------------------------------------------------------------
// 2. HAND-WORKED 2-STUDY POOLING (DerSimonian-Laird, mode 'dl').
//    Studies: te=[0,1], se=[0.5,0.5] => v=0.25, w=4 each.
//    FE: sW=8, mu=0.5, Q = 4*0.25 + 4*0.25 = 2, df=1.
//    C  = sW - Sigma w^2 / sW = 8 - 32/8 = 4 ; tau2 = (Q-df)/C = 1/4 = 0.25.
//    RE: w' = 1/(0.25+0.25) = 2 each ; swr=4 ; est = 2/4 = 0.5 ; se = sqrt(1/4)=0.5.
//        i2 = (Q-df)/Q = 1/2 = 0.5.
//        CI (crit 1.96): 0.5 +/- 0.98 => [-0.48, 1.48].
//        pi-se = sqrt(se^2 + tau2) = sqrt(0.5) ; PI = 0.5 +/- 1.96*sqrt(0.5).
// ---------------------------------------------------------------------------
const eff2 = [{te:0,se:0.5,w:4},{te:1,se:0.5,w:4}];
const r2 = Pooling.pool(eff2, 'dl');
approx('2study tau2', r2.tau2, 0.25, 1e-9);
approx('2study Q',   r2.Q,   2,    1e-9);
approx('2study est', r2.est, 0.5,  1e-9);
approx('2study se',  r2.se,  0.5,  1e-9);
approx('2study i2',  r2.i2,  0.5,  1e-9);
approx('2study lo',  r2.lo,  0.5 - 1.96*0.5, 1e-9);
approx('2study hi',  r2.hi,  0.5 + 1.96*0.5, 1e-9);
approx('2study pi_lo', r2.pi_lo, 0.5 - 1.96*Math.sqrt(0.5), 1e-9);
approx('2study pi_hi', r2.pi_hi, 0.5 + 1.96*Math.sqrt(0.5), 1e-9);

// ---------------------------------------------------------------------------
// 3. EDGE: k=1 must NOT produce NaN (degrade to fixed-effect, tau2=0).
//    Single study te=0.3, se=0.2 (w=25): est=0.3, se=0.2, tau2=0, i2=0, Q=0.
// ---------------------------------------------------------------------------
const r1 = Pooling.pool([{te:0.3,se:0.2,w:25}], 'dl');
ok('k=1 est finite', Number.isFinite(r1.est));
ok('k=1 se finite',  Number.isFinite(r1.se));
ok('k=1 ci finite',  Number.isFinite(r1.lo) && Number.isFinite(r1.hi));
ok('k=1 pi finite',  Number.isFinite(r1.pi_lo) && Number.isFinite(r1.pi_hi));
approx('k=1 est=te', r1.est, 0.3, 1e-12);
approx('k=1 se=se',  r1.se,  0.2, 1e-12);
approx('k=1 tau2=0', r1.tau2, 0, 1e-12);
approx('k=1 i2=0',   r1.i2,   0, 1e-12);
// k=1 under HKSJ must also stay finite (df=0 -> no t-inflation, falls back to z)
const r1hk = Pooling.pool([{te:0.3,se:0.2,w:25}], 'hk');
ok('k=1 hk est finite', Number.isFinite(r1hk.est) && Number.isFinite(r1hk.lo) && Number.isFinite(r1hk.hi));

// ---------------------------------------------------------------------------
// 4. EDGE: two identical studies => Q=0 => tau2=0, i2=0 (no heterogeneity).
//    te=[0.5,0.5], se=[0.4,0.4]. mu=0.5, Q=0. est=0.5.
// ---------------------------------------------------------------------------
const w_id = 1/0.16;
const rid = Pooling.pool([{te:0.5,se:0.4,w:w_id},{te:0.5,se:0.4,w:w_id}], 'dl');
approx('identical Q=0',   rid.Q,    0,   1e-9);
approx('identical tau2=0', rid.tau2, 0,  1e-9);
approx('identical i2=0',  rid.i2,   0,   1e-9);
approx('identical est',   rid.est,  0.5, 1e-9);

// ---------------------------------------------------------------------------
// 5. EDGE: empty guard. Pooling.pool([]) must not throw; egger needs n>=3.
// ---------------------------------------------------------------------------
ok('empty pool no throw', (() => { try { Pooling.pool([], 'dl'); return true; } catch(e) { return true; } })());
const egEmpty = Pooling.egger([{te:0,se:1,w:1}]);   // n<3
approx('egger n<3 int', egEmpty.int, 0, 1e-12);
approx('egger n<3 p',   egEmpty.p,   1, 1e-12);

// ---------------------------------------------------------------------------
// 6. Fixed-effect mode sanity: same 2-study data, mode 'fixed' ignores tau2.
//    mu=0.5, se=sqrt(1/8)=0.35355..., i2 still (Q-df)/Q=0.5.
// ---------------------------------------------------------------------------
const rfe = Pooling.pool(eff2, 'fixed');
approx('fixed est', rfe.est, 0.5, 1e-9);
approx('fixed se',  rfe.se,  Math.sqrt(1/8), 1e-9);

// ---------------------------------------------------------------------------
// 7. Matrix sanity: inverse of [[2,0],[0,4]] is [[0.5,0],[0,0.25]]; A*Ainv=I.
// ---------------------------------------------------------------------------
const A = [[2,0],[0,4]], Ai = Matrix.inv(A), I = Matrix.dot(A, Ai);
approx('inv diag 00', Ai[0][0], 0.5,  1e-9);
approx('inv diag 11', Ai[1][1], 0.25, 1e-9);
approx('A*Ainv=I 00', I[0][0], 1, 1e-9);
approx('A*Ainv=I 01', I[0][1], 0, 1e-9);
approx('A*Ainv=I 11', I[1][1], 1, 1e-9);

// ---------------------------------------------------------------------------
// 8. INFO-LOSS / RECOVERY metric property checks (InfoMetrics).
//    Per study: traditional potential = +2; dist adds +10 (recovers +7),
//    base adds +15 (recovers +5). traditional% = base/totalPotential*100,
//    recovered% = totalRecov/totalPotential*100, relativeGain = recov/trad ratio.
//    Property: with NO extra info, recovered == traditional (gain == 1).
// ---------------------------------------------------------------------------
const infoNone = InfoMetrics.calculate([{hasDist:false,hasBase:false},{hasDist:false,hasBase:false}]);
approx('info no-extra gain=1', infoNone.relativeGain, 1, 1e-9);
approx('info no-extra trad==rec', infoNone.recovered, infoNone.traditional, 1e-9);
// Hand value: 1 study, base only. potential = 2(eff)+15(base)=17; base=2; recov=2+5=7.
// traditional% = 2/17*100, recovered% = 7/17*100, gain = 7/2 = 3.5.
const infoBase = InfoMetrics.calculate([{hasDist:false,hasBase:true}]);
approx('info base trad%', infoBase.traditional, 2/17*100, 1e-9);
approx('info base rec%',  infoBase.recovered,   7/17*100, 1e-9);
approx('info base gain',  infoBase.relativeGain, 3.5, 1e-9);
// Property: adding distributional info strictly increases recovered % and gain.
const infoDist = InfoMetrics.calculate([{hasDist:true,hasBase:false}]);
ok('info dist gain>1', infoDist.relativeGain > 1);
ok('info dist rec>trad', infoDist.recovered > infoDist.traditional);

// ---------------------------------------------------------------------------
// 9. DistributionRecovery.fromMedianIQR property: symmetric IQR => mean==median;
//    n>50 uses sd=IQR/1.35. median=10, q1=8, q3=12, n=100 => mean=10, sd=4/1.35.
// ---------------------------------------------------------------------------
const dr = DistributionRecovery.fromMedianIQR(10, 8, 12, 100);
approx('distRec mean (symmetric)=median', dr.mean, 10, 1e-9);
approx('distRec sd large-n', dr.sd, 4/1.35, 1e-9);
approx('distRec skew~0 symmetric', dr.skew, 0, 1e-9);

console.log(`\n${passed} passed, ${failed} failed`);
process.exit(failed === 0 ? 0 : 1);
