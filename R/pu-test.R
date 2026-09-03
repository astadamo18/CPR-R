# Phillips-Ouliaris-type PU cointegration test for a cointegrating
# polynomial regression. Port of PU_test.m.
#
# Unlike ct_test() (which needs FM-OLS residuals and Omega_udotv1 from a
# fitted cpr() object), PU_test operates directly on the raw y/x series: it
# runs its own plain-OLS regression to get residuals for the denominator,
# and a separate local VAR(1)-with-deterministics on the stacked [y, x]
# system for the long-run variance in the numerator. So it is not simply
# "the same computation with different critical values" -- it is a
# genuinely different test statistic, reusing the same kernel/bandwidth/
# long-run-variance machinery as cpr() and ct_test().
#
# Rewritten the same way ct_test() was: an S3 generic dispatching either on
# the raw ingredients (pu_test.default()) or directly on a fitted cpr()
# object (pu_test.cpr(), which pulls y/x/d/m/p off the fit -- possible
# because cpr() now stores the resolved raw y/x on the returned object),
# plus the full critical-value grid (all 48 (d, m, p) combinations, from
# every PUcritval/*.mat file in the original toolbox) so the matching table
# is picked automatically instead of the caller managing one by hand.

.pu_critval_table <- list(
  "d-1_m1_p1" = c(`1%` = 0.5844927244472804, `2.5%` = 0.8326050258687271, `5%` = 1.1653479660402613, `10%` = 1.7664537722022524, `50%` = 7.229156664205592, `90%` = 20.50395670086267, `95%` = 26.11959425621969, `97.5%` = 31.652640770302234, `99%` = 39.040326164205915),
  "d-1_m1_p2" = c(`1%` = 0.8033262394763565, `2.5%` = 1.1260638676048425, `5%` = 1.5464587084769665, `10%` = 2.288931388746666, `50%` = 8.957651912745531, `90%` = 24.504372015553876, `95%` = 30.87544646350799, `97.5%` = 36.98618345968528, `99%` = 44.86187550260252),
  "d-1_m1_p3" = c(`1%` = 0.9509483085706806, `2.5%` = 1.3361580117703347, `5%` = 1.8349966159931852, `10%` = 2.6903518507250346, `50%` = 10.072428787191695, `90%` = 27.079380307902603, `95%` = 33.85960756217187, `97.5%` = 40.5754747038345, `99%` = 48.977050701870205),
  "d-1_m1_p4" = c(`1%` = 1.0791400534619768, `2.5%` = 1.5078584565749893, `5%` = 2.0546263038447936, `10%` = 3.0064025804123586, `50%` = 10.866285030127983, `90%` = 28.880146540379467, `95%` = 36.18056714920501, `97.5%` = 43.20610925679325, `99%` = 51.80405022442759),
  "d-1_m2_p1" = c(`1%` = 1.0623060138836977, `2.5%` = 1.604802973069037, `5%` = 2.298469200961652, `10%` = 3.443374188273747, `50%` = 11.220315393333765, `90%` = 27.249937258290466, `95%` = 33.512196775443094, `97.5%` = 39.815550162107776, `99%` = 47.73408153656163),
  "d-1_m2_p2" = c(`1%` = 1.3403722707091044, `2.5%` = 1.9742577671493147, `5%` = 2.790217784529129, `10%` = 4.147533727568113, `50%` = 13.295171372256114, `90%` = 31.169902335863316, `95%` = 38.03631773163201, `97.5%` = 44.871608985117746, `99%` = 53.22644064641601),
  "d-1_m2_p3" = c(`1%` = 1.5754925379488463, `2.5%` = 2.292837533855396, `5%` = 3.1753387831550737, `10%` = 4.641231467354918, `50%` = 14.595847073661686, `90%` = 33.822578985719915, `95%` = 41.22162887760914, `97.5%` = 48.374018070980576, `99%` = 57.010242346940174),
  "d-1_m2_p4" = c(`1%` = 1.7339019537180467, `2.5%` = 2.4980776557640043, `5%` = 3.4702728726612744, `10%` = 5.015880614651444, `50%` = 15.525202746908105, `90%` = 35.82401359796987, `95%` = 43.64275853280451, `97.5%` = 51.04641337468866, `99%` = 59.95098559557164),
  "d-1_m3_p1" = c(`1%` = 1.8991407805397928, `2.5%` = 2.808694367518194, `5%` = 3.907300830432876, `10%` = 5.478813895509526, `50%` = 15.26634112242204, `90%` = 33.35443914625714, `95%` = 40.32132802914055, `97.5%` = 47.03898359511122, `99%` = 55.694772168474834),
  "d-1_m3_p2" = c(`1%` = 2.2619069223069648, `2.5%` = 3.340653026808762, `5%` = 4.574431863623435, `10%` = 6.389082073402556, `50%` = 17.54052675298741, `90%` = 37.36044783816451, `95%` = 44.74090616549473, `97.5%` = 51.650311144121645, `99%` = 60.93165835770328),
  "d-1_m3_p3" = c(`1%` = 2.5319784902655194, `2.5%` = 3.7200837332967165, `5%` = 5.026235406282889, `10%` = 6.9524608983601786, `50%` = 18.98653859614306, `90%` = 40.09453519566864, `95%` = 47.82026950783646, `97.5%` = 55.22595566272133, `99%` = 64.87981232474867),
  "d-1_m3_p4" = c(`1%` = 2.767664491402295, `2.5%` = 4.008294942787907, `5%` = 5.360045293555596, `10%` = 7.394351315861529, `50%` = 20.03525724524044, `90%` = 42.10770109640377, `95%` = 50.30663049525595, `97.5%` = 57.86136852885052, `99%` = 67.73637863640467),
  "d-1_m4_p1" = c(`1%` = 3.0652723287569517, `2.5%` = 4.319259755755328, `5%` = 5.72714572677334, `10%` = 7.810661396181103, `50%` = 19.19635810704214, `90%` = 39.05500997115648, `95%` = 46.59911694419849, `97.5%` = 53.625813747299006, `99%` = 62.52839675105327),
  "d-1_m4_p2" = c(`1%` = 3.501756133834823, `2.5%` = 4.98686143988409, `5%` = 6.642714535534358, `10%` = 8.971854289040536, `50%` = 21.692853248621436, `90%` = 43.114368515718354, `95%` = 51.00964757949137, `97.5%` = 58.70393517466559, `99%` = 68.03975966284908),
  "d-1_m4_p3" = c(`1%` = 3.816501668714709, `2.5%` = 5.445827697140178, `5%` = 7.179909025056194, `10%` = 9.635360191763002, `50%` = 23.320401628313835, `90%` = 45.95773624130684, `95%` = 54.19440087411578, `97.5%` = 61.996261607723284, `99%` = 71.79955496204273),
  "d-1_m4_p4" = c(`1%` = 4.1028524328481, `2.5%` = 5.75033909435192, `5%` = 7.585643673111426, `10%` = 10.147394385854293, `50%` = 24.4741403648446, `90%` = 48.03325457581432, `95%` = 56.629539607351234, `97.5%` = 64.64508330900254, `99%` = 74.75732872297571),
  "d0_m1_p1" = c(`1%` = 1.8824999375854483, `2.5%` = 2.475231916943531, `5%` = 3.179045332564682, `10%` = 4.311631387546473, `50%` = 11.969602287611671, `90%` = 27.77892816145555, `95%` = 34.097033795530244, `97.5%` = 40.14736646973793, `99%` = 48.03478123456876),
  "d0_m1_p2" = c(`1%` = 2.0262593733019623, `2.5%` = 2.6460146295490197, `5%` = 3.4310662947470956, `10%` = 4.6943570569576085, `50%` = 13.436969556179534, `90%` = 30.982539242290752, `95%` = 37.87475517642149, `97.5%` = 44.432484707331454, `99%` = 52.94784857075556),
  "d0_m1_p3" = c(`1%` = 2.0990719611340376, `2.5%` = 2.7548982591790656, `5%` = 3.5723525409246784, `10%` = 4.8906700938804155, `50%` = 14.279660307794828, `90%` = 33.153798502437304, `95%` = 40.436966766820916, `97.5%` = 47.47759570076783, `99%` = 56.12298638144296),
  "d0_m1_p4" = c(`1%` = 2.144165488949347, `2.5%` = 2.8180854459262608, `5%` = 3.6671994761157083, `10%` = 5.025606032979009, `50%` = 14.849545242069196, `90%` = 34.71771710014963, `95%` = 42.42543925311394, `97.5%` = 49.79570721133383, `99%` = 58.788823187382626),
  "d0_m2_p1" = c(`1%` = 2.581572968486881, `2.5%` = 3.498372591124679, `5%` = 4.58439576323307, `10%` = 6.15753419700135, `50%` = 15.765584118523488, `90%` = 33.84948428586193, `95%` = 40.903922096590584, `97.5%` = 47.677081928780524, `99%` = 55.79960962394327),
  "d0_m2_p2" = c(`1%` = 2.786751917495254, `2.5%` = 3.792536903918538, `5%` = 5.018915807414055, `10%` = 6.760052634986468, `50%` = 17.61292281600591, `90%` = 37.27683093053823, `95%` = 44.791601010864454, `97.5%` = 51.84105938287113, `99%` = 60.60539570189891),
  "d0_m2_p3" = c(`1%` = 2.8737447784499337, `2.5%` = 3.9536033992546717, `5%` = 5.240956890194061, `10%` = 7.0984347222346456, `50%` = 18.705422676972546, `90%` = 39.67625680572178, `95%` = 47.70207827233582, `97.5%` = 54.92065818431272, `99%` = 63.828305873786796),
  "d0_m2_p4" = c(`1%` = 2.946667262401094, `2.5%` = 4.043021619615831, `5%` = 5.385205400699162, `10%` = 7.309631856277119, `50%` = 19.454339634960647, `90%` = 41.4570995662834, `95%` = 49.662834943641336, `97.5%` = 57.29633787736984, `99%` = 66.69610310585738),
  "d0_m3_p1" = c(`1%` = 3.6827861428337747, `2.5%` = 4.914305169917288, `5%` = 6.292296653662086, `10%` = 8.262153127559206, `50%` = 19.679411518620057, `90%` = 39.52474073101097, `95%` = 46.892009462823076, `97.5%` = 53.84270194677741, `99%` = 63.170581742789885),
  "d0_m3_p2" = c(`1%` = 3.9874386348749806, `2.5%` = 5.341249876551831, `5%` = 6.881057898405195, `10%` = 9.105518686225956, `50%` = 21.713717518777894, `90%` = 43.08627009862985, `95%` = 50.91644221882318, `97.5%` = 58.109997712873074, `99%` = 67.8341506076325),
  "d0_m3_p3" = c(`1%` = 4.143774103437503, `2.5%` = 5.575765670470416, `5%` = 7.19042203361139, `10%` = 9.554223064847058, `50%` = 23.013064932115896, `90%` = 45.569984257762904, `95%` = 53.63964764336589, `97.5%` = 61.309760644011504, `99%` = 71.14782413507763),
  "d0_m3_p4" = c(`1%` = 4.270276041528504, `2.5%` = 5.728784438116732, `5%` = 7.392169986806199, `10%` = 9.855537081136518, `50%` = 23.920973504147177, `90%` = 47.46618663770211, `95%` = 55.869175990667145, `97.5%` = 63.83567098991168, `99%` = 73.93766985402081),
  "d0_m4_p1" = c(`1%` = 4.918606439308174, `2.5%` = 6.493747463812977, `5%` = 8.224468972595169, `10%` = 10.624793479858264, `50%` = 23.57289552309916, `90%` = 45.16688632694942, `95%` = 52.994986157968135, `97.5%` = 60.256991018052474, `99%` = 69.95229166470926),
  "d0_m4_p2" = c(`1%` = 5.318790623828896, `2.5%` = 7.123715071012359, `5%` = 9.06055699215621, `10%` = 11.687640151011964, `50%` = 25.84050790097934, `90%` = 48.73910935877864, `95%` = 56.92444655137675, `97.5%` = 64.94504171418026, `99%` = 74.38936336194989),
  "d0_m4_p3" = c(`1%` = 5.557877289152626, `2.5%` = 7.4697108163994805, `5%` = 9.498807880621863, `10%` = 12.269472661809557, `50%` = 27.331753242934955, `90%` = 51.37814216811716, `95%` = 59.99154925283472, `97.5%` = 68.03690746569595, `99%` = 78.35666278593254),
  "d0_m4_p4" = c(`1%` = 5.683748615855976, `2.5%` = 7.654529541292879, `5%` = 9.797859742014104, `10%` = 12.655998745503647, `50%` = 28.357286921565457, `90%` = 53.370219575676295, `95%` = 62.17429712826514, `97.5%` = 70.4036910113825, `99%` = 81.1152845336433),
  "d1_m1_p1" = c(`1%` = 5.516346391948592, `2.5%` = 6.802259607232868, `5%` = 8.192782178891173, `10%` = 10.151519049624655, `50%` = 21.433180726548063, `90%` = 41.13153495107404, `95%` = 48.563955131764246, `97.5%` = 55.32907458832227, `99%` = 64.46303899796929),
  "d1_m1_p2" = c(`1%` = 6.280123912958905, `2.5%` = 7.7520727699108525, `5%` = 9.336522932742819, `10%` = 11.582632288651025, `50%` = 24.081925409651433, `90%` = 45.237292689766434, `95%` = 52.95181044105902, `97.5%` = 60.32064282500815, `99%` = 69.4705880811974),
  "d1_m1_p3" = c(`1%` = 6.557617881447265, `2.5%` = 8.11817770928737, `5%` = 9.81204692721518, `10%` = 12.259223933629926, `50%` = 25.55267182099244, `90%` = 47.92519672341596, `95%` = 55.92644145664932, `97.5%` = 63.59797606346454, `99%` = 73.3073809641591),
  "d1_m1_p4" = c(`1%` = 6.73535029548124, `2.5%` = 8.383276894066505, `5%` = 10.15006833389631, `10%` = 12.710214161862348, `50%` = 26.62096573500015, `90%` = 50.014276841686616, `95%` = 58.27906900065073, `97.5%` = 66.04700040621623, `99%` = 76.12637297612757),
  "d1_m2_p1" = c(`1%` = 6.544413369813812, `2.5%` = 8.19677307289401, `5%` = 9.837650650828724, `10%` = 12.157881305021824, `50%` = 25.044553790718265, `90%` = 46.396681251828866, `95%` = 54.26958661293958, `97.5%` = 61.80348227049035, `99%` = 71.49864501167421),
  "d1_m2_p2" = c(`1%` = 7.463571110493297, `2.5%` = 9.246260154771756, `5%` = 11.16004916706096, `10%` = 13.75057175290935, `50%` = 27.800208461122587, `90%` = 50.424206943994356, `95%` = 58.667259138665166, `97.5%` = 66.38472191939523, `99%` = 76.36745975607384),
  "d1_m2_p3" = c(`1%` = 7.798597637598878, `2.5%` = 9.701226986741938, `5%` = 11.689941803262776, `10%` = 14.49388515644401, `50%` = 29.45148356077844, `90%` = 53.36317956013909, `95%` = 61.84378630130338, `97.5%` = 69.85482288291209, `99%` = 79.93884142771154),
  "d1_m2_p4" = c(`1%` = 8.050589471743907, `2.5%` = 10.05505396700691, `5%` = 12.09660543752007, `10%` = 15.036785701394566, `50%` = 30.652126242246258, `90%` = 55.39044400235098, `95%` = 64.384124299302, `97.5%` = 72.57864449108914, `99%` = 82.62721524315761),
  "d1_m3_p1" = c(`1%` = 7.855627934266504, `2.5%` = 9.783904340469288, `5%` = 11.77273167290594, `10%` = 14.546696565246798, `50%` = 28.801249497923845, `90%` = 51.71619380407528, `95%` = 59.874319272960214, `97.5%` = 67.58522396435923, `99%` = 77.26377510203442),
  "d1_m3_p2" = c(`1%` = 8.85610866898533, `2.5%` = 11.018747060631028, `5%` = 13.180904742060909, `10%` = 16.231094910331333, `50%` = 31.66901812624927, `90%` = 55.63683951585686, `95%` = 64.25646157332216, `97.5%` = 72.25362329179605, `99%` = 82.47735109559927),
  "d1_m3_p3" = c(`1%` = 9.252877905895152, `2.5%` = 11.567097254474476, `5%` = 13.919405105943458, `10%` = 17.164344298855244, `50%` = 33.419960581668306, `90%` = 58.44192338140805, `95%` = 67.32719971830115, `97.5%` = 75.53420166139988, `99%` = 86.6034616254694),
  "d1_m3_p4" = c(`1%` = 9.535180590996537, `2.5%` = 11.924877652910975, `5%` = 14.422529876950458, `10%` = 17.796527112127162, `50%` = 34.70566731237989, `90%` = 60.640849963945705, `95%` = 69.76756441450567, `97.5%` = 78.15902592159853, `99%` = 89.63867219597928),
  "d1_m4_p1" = c(`1%` = 9.191454344614808, `2.5%` = 11.401089942697999, `5%` = 13.700442596233916, `10%` = 16.83826304437877, `50%` = 32.48192118900012, `90%` = 56.804619813094355, `95%` = 65.3239939487248, `97.5%` = 73.24460201755508, `99%` = 83.55020545022833),
  "d1_m4_p2" = c(`1%` = 10.351890664592599, `2.5%` = 12.73932314744644, `5%` = 15.225335370035264, `10%` = 18.644766838575524, `50%` = 35.348672367556944, `90%` = 60.89185700900691, `95%` = 69.89283607168073, `97.5%` = 78.24591558414212, `99%` = 89.02332114503224),
  "d1_m4_p3" = c(`1%` = 10.857732894240154, `2.5%` = 13.40665237599865, `5%` = 16.07109272355315, `10%` = 19.671529020989407, `50%` = 37.27519843969966, `90%` = 63.91134227365861, `95%` = 73.03919331012547, `97.5%` = 81.61911030502118, `99%` = 92.44382980280136),
  "d1_m4_p4" = c(`1%` = 11.213422029239261, `2.5%` = 13.86342955537533, `5%` = 16.60477351604415, `10%` = 20.373831562325687, `50%` = 38.65323472139553, `90%` = 66.15056980338596, `95%` = 75.55621324702476, `97.5%` = 84.44006900777497, `99%` = 95.65022834463875)
)

#' @keywords internal
pu_critval <- function(d, m, p) {
  key <- paste0("d", d, "_m", m, "_p", p)
  tab <- .pu_critval_table[[key]]
  if (is.null(tab)) {
    stop("No PU critical value table for d=", d, ", m=", m, ", p=", p,
         ". The bundled grid covers d in {-1, 0, 1}, m in {1, 2, 3, 4}, ",
         "p in {1, 2, 3, 4} (the full range the original CTPUTests toolbox ",
         "itself tabulates); this combination falls outside it.", call. = FALSE)
  }
  tab
}

#' Phillips-Ouliaris-type PU cointegration test for a cointegrating polynomial regression
#'
#' Tests the null hypothesis of *no* cointegration against the alternative
#' of cointegration (the reverse null from [ct_test()]), using the
#' Phillips-Ouliaris/Wagner-type ratio statistic of Wagner and co-authors
#' for cointegrating polynomial regressions. Port of `PU_test.m`.
#'
#' An S3 generic, mirroring [ct_test()]: call it either on the raw
#' ingredients (`pu_test(y, x, d, m, orders, kernel, bandwidth, alpha)`,
#' the [pu_test.default()] method) or directly on a fitted [cpr()] object
#' (`pu_test(fit, alpha = ...)`, the [pu_test.cpr()] method), which pulls
#' `y`, `x`, `m`, `orders`, `kernel`, `bandwidth`, and `d` off the fit for
#' you. Either way, the matching critical value table for the resulting
#' `(d, m, p)` is picked automatically from the full bundled grid (`d` in
#' `{-1, 0, 1}`, `m` and `p` in `{1, 2, 3, 4}`).
#'
#' Unlike `ct_test()`, this does not use the fit's FM-OLS residuals: it
#' runs its own plain-OLS regression internally for the denominator, and a
#' separate local VAR(1)-with-deterministics on the stacked `[y, x]` system
#' for the long-run variance in the numerator -- so it works off any `cpr`
#' fit regardless of `estimator`.
#'
#' The returned object has a `print()` method showing the test statistic,
#' critical values, decisions, and hypotheses (like [ct_test()]'s).
#'
#' @param y Either the dependent variable (length `T`, dispatching to
#'   [pu_test.default()]) or a fitted [cpr()] object (dispatching to
#'   [pu_test.cpr()]).
#' @param ... Passed on to the method.
#'
#' @return An object of class `"pu_test"`: a list with `statistic`,
#'   `alpha`, `critval` (critical value per `alpha`), `reject` (logical per
#'   `alpha`; `TRUE` means reject the null of no cointegration, i.e.
#'   evidence *for* cointegration), and `d`/`m`/`p` (the critical value
#'   table used).
#' @export
pu_test <- function(y, ...) {
  UseMethod("pu_test")
}

#' @describeIn pu_test Default method: supply the series and specification
#'   directly, as in `PU_test.m`.
#' @param x Integrated regressors, `T x m`. If `m > 1`, only the *last*
#'   column gets the polynomial powers requested via `orders`; the
#'   remaining `m - 1` columns enter linearly (power 1).
#' @param d Deterministic specification: `-1` (none), `0` (intercept), `1`
#'   (intercept + trend).
#' @param m Number of integrated regressors (`ncol(x)`).
#' @param orders Powers of the last column of `x` to include; see
#'   [gen_var_poly_terms()]. If `orders` is exactly `1` (linear only), no
#'   critical values exist in the original toolbox either, and
#'   `critval`/`reject` come back as `NA`.
#' @param kernel Kernel function, see [lr_weights()].
#' @param bandwidth Bandwidth selection: `"And91"`, `"AM92"`, `"NW"`, or a
#'   fixed numeric value.
#' @param alpha Significance levels to test at; must be a subset of
#'   `c(0.1, 0.05, 0.01)`.
#' @export
pu_test.default <- function(y, x, d, m, orders, kernel, bandwidth, alpha = c(0.1, 0.05, 0.01), ...) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  Tn <- nrow(y)

  const <- rep(1, Tn)
  trend <- seq_len(Tn)
  deter <- switch(as.character(d),
    "-1" = matrix(numeric(0), Tn, 0),
    "0"  = matrix(const, ncol = 1, dimnames = list(NULL, "const")),
    "1"  = cbind(const = const, trend = trend),
    stop("`d` must be -1, 0, or 1.", call. = FALSE)
  )

  z <- cbind(y, x)

  linear_only <- is.numeric(orders) && length(orders) == 1 && orders == 1
  xpower <- x[, m, drop = FALSE]

  if (linear_only) {
    powerreg <- NULL
  } else {
    powerreg <- gen_var_poly_terms(xpower, orders, stochastic = FALSE)$X
  }

  if (linear_only) {
    regmat4uhat <- cbind(deter, x)
  } else {
    x_rest <- if (m > 1) x[, seq_len(m - 1), drop = FALSE] else NULL
    regmat4uhat <- cbind(deter, x_rest, powerreg)
  }

  coeff4uhat <- solve(crossprod(regmat4uhat), crossprod(regmat4uhat, y))
  uhat <- y - regmat4uhat %*% coeff4uhat

  # Local VAR(1)-with-deterministics on the stacked [y, x] system:
  depvar <- z[2:Tn, , drop = FALSE]
  lagz <- lag_matrix(z, 1)[2:Tn, , drop = FALSE]
  indepvar <- cbind(deter[2:Tn, , drop = FALSE], lagz)

  varcoeff <- solve(crossprod(indepvar), crossprod(indepvar, depvar))
  varresid <- depvar - indepvar %*% varcoeff

  Lr <- estimate_lr_var(varresid, kernel, bandwidth, demean = FALSE)$Omega
  omega_udotv <- as.numeric(Lr[1, 1] - Lr[1, -1, drop = FALSE] %*% solve(Lr[-1, -1, drop = FALSE], Lr[-1, 1]))

  statistic <- omega_udotv / (Tn^(-2) * sum(uhat^2))

  if (linear_only) {
    return(structure(
      list(statistic = statistic, alpha = alpha,
           critval = rep(NA_real_, length(alpha)), reject = rep(NA, length(alpha)),
           d = d, m = m, p = orders),
      class = "pu_test"
    ))
  }

  crit_table <- pu_critval(d, m, orders)
  pct_label <- c(`0.1` = "90%", `0.05` = "95%", `0.01` = "99%")
  lab <- pct_label[as.character(alpha)]
  if (any(is.na(lab))) {
    stop("`alpha` must be a subset of c(0.1, 0.05, 0.01).", call. = FALSE)
  }
  critval <- as.numeric(crit_table[lab])

  structure(
    list(statistic = statistic, alpha = alpha, critval = critval,
         reject = statistic > critval, d = d, m = m, p = orders),
    class = "pu_test"
  )
}

#' @describeIn pu_test `cpr` method: `y`, `x`, `kernel`, and `bandwidth` are
#'   read straight off the fit (`x$y`, `x$x`, `x$kernel`, `x$bandwidth`);
#'   `d` defaults to `NULL`, which infers it from the fit's own `deter`
#'   (the same way [ct_test.cpr()] does). Only supports a fit with a single
#'   integrated regressor whose powers are the sequential `1:p` (i.e. the
#'   default `orders = p` form) -- the structural case `PU_test.m` itself
#'   supports via a single `orders` scalar; for anything else (multiple
#'   regressors, non-sequential powers), call [pu_test.default()] directly
#'   with the right `d`/`m`/`orders`.
#' @export
pu_test.cpr <- function(y, d = NULL, alpha = c(0.1, 0.05, 0.01), ...) {
  fit <- y
  y_raw <- fit$y
  x_raw <- fit$x
  if (is.null(y_raw) || is.null(x_raw)) {
    stop("This cpr fit does not carry the raw y/x pu_test() needs (re-fit with a ",
         "current version of cpr(), which stores them on the returned object).",
         call. = FALSE)
  }

  m <- ncol(x_raw)
  powers1 <- fit$fit$powers[[1]]
  if (m != 1 || !identical(powers1, seq_len(max(powers1)))) {
    stop("pu_test() on a cpr fit only supports a single integrated regressor with ",
         "sequential powers 1:p (this fit has ", m, " regressor(s)",
         if (m == 1) paste0(" with powers ", paste(powers1, collapse = ", ")) else "",
         "). Call pu_test.default(y, x, d, m, orders, kernel, bandwidth) directly ",
         "with the right specification instead.", call. = FALSE)
  }
  p <- max(powers1)

  if (is.null(d)) d <- infer_cpr_d(fit)

  pu_test.default(y_raw, x_raw, d = d, m = m, orders = p,
                   kernel = fit$kernel, bandwidth = fit$bandwidth, alpha = alpha)
}

#' @export
print.pu_test <- function(x, digits = 4, ...) {
  cat("PU test for cointegration\n\n")
  cat("H0: no cointegration        H1: cointegration\n\n")
  cat("d =", x$d, " m =", x$m, " p =", x$p, "\n")
  cat("Test statistic:", format(round(x$statistic, digits)), "\n\n")

  alpha_pct <- paste0(format(100 * x$alpha, trim = TRUE), "%")
  cv <- x$critval
  names(cv) <- alpha_pct
  decision <- ifelse(is.na(x$reject), "n/a", ifelse(x$reject, "reject", "do not reject"))
  names(decision) <- alpha_pct

  cat("Critical values:\n")
  print(round(cv, digits))
  cat("\nDecision:\n")
  print(decision, quote = FALSE)
  invisible(x)
}
