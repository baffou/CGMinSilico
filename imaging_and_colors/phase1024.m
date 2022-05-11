function phasemap=phase1024(npx)

phasemap = [...
1 1 1
1 0.99951 0.99951
0.99988 0.99708 0.99705
0.99884 0.99466 0.9946
0.9978 0.99225 0.99216
0.99676 0.98986 0.98971
0.99571 0.98748 0.98727
0.99467 0.98512 0.98483
0.99363 0.98277 0.9824
0.99259 0.98044 0.97997
0.99155 0.97812 0.97754
0.99051 0.97582 0.97511
0.98947 0.97353 0.97269
0.98843 0.97125 0.97027
0.98739 0.96899 0.96785
0.98634 0.96675 0.96544
0.9853 0.96451 0.96303
0.98426 0.9623 0.96062
0.98322 0.96009 0.95821
0.98218 0.9579 0.95581
0.98114 0.95573 0.95341
0.9801 0.95357 0.95101
0.97906 0.95142 0.94862
0.97801 0.94928 0.94623
0.97697 0.94716 0.94384
0.97593 0.94506 0.94146
0.97489 0.94297 0.93907
0.97385 0.94089 0.93669
0.97281 0.93882 0.93432
0.97177 0.93677 0.93195
0.97073 0.93474 0.92958
0.96969 0.93271 0.92721
0.96864 0.9307 0.92484
0.9676 0.92871 0.92248
0.96656 0.92672 0.92013
0.96552 0.92475 0.91777
0.96448 0.9228 0.91542
0.96344 0.92086 0.91307
0.9624 0.91893 0.91072
0.96136 0.91701 0.90838
0.96032 0.91511 0.90604
0.95927 0.91322 0.9037
0.95823 0.91135 0.90136
0.95719 0.90948 0.89903
0.95615 0.90764 0.8967
0.95511 0.9058 0.89438
0.95407 0.90398 0.89206
0.95303 0.90217 0.88974
0.95199 0.90037 0.88742
0.95095 0.89859 0.8851
0.9499 0.89682 0.88279
0.94886 0.89506 0.88049
0.94782 0.89331 0.87818
0.94678 0.89158 0.87588
0.94574 0.88986 0.87358
0.9447 0.88816 0.87128
0.94366 0.88646 0.86899
0.94262 0.88478 0.8667
0.94158 0.88311 0.86441
0.94053 0.88146 0.86213
0.93949 0.87982 0.85984
0.93845 0.87819 0.85757
0.93741 0.87657 0.85529
0.93637 0.87496 0.85302
0.93533 0.87337 0.85075
0.93429 0.87179 0.84848
0.93325 0.87022 0.84622
0.93221 0.86866 0.84396
0.93116 0.86712 0.8417
0.93012 0.86559 0.83944
0.92908 0.86407 0.83719
0.92804 0.86256 0.83494
0.927 0.86107 0.83269
0.92596 0.85959 0.83045
0.92492 0.85812 0.82821
0.92388 0.85666 0.82597
0.92284 0.85521 0.82374
0.92179 0.85378 0.82151
0.92075 0.85236 0.81928
0.91971 0.85095 0.81705
0.91867 0.84955 0.81483
0.91763 0.84816 0.81261
0.91659 0.84678 0.81039
0.91555 0.84542 0.80818
0.91451 0.84407 0.80597
0.91347 0.84273 0.80376
0.91242 0.8414 0.80156
0.91138 0.84008 0.79936
0.91034 0.83878 0.79716
0.9093 0.83749 0.79496
0.90826 0.8362 0.79277
0.90722 0.83493 0.79058
0.90618 0.83367 0.78839
0.90514 0.83243 0.7862
0.9041 0.83119 0.78402
0.90305 0.82997 0.78184
0.90201 0.82875 0.77967
0.90097 0.82755 0.7775
0.89993 0.82636 0.77533
0.89889 0.82518 0.77316
0.89785 0.82401 0.771
0.89681 0.82285 0.76884
0.89577 0.8217 0.76668
0.89472 0.82057 0.76452
0.89368 0.81944 0.76237
0.89264 0.81833 0.76022
0.8916 0.81722 0.75808
0.89056 0.81613 0.75593
0.88952 0.81505 0.75379
0.88848 0.81398 0.75166
0.88744 0.81292 0.74952
0.8864 0.81187 0.74739
0.88535 0.81083 0.74526
0.88431 0.8098 0.74314
0.88338 0.80811 0.74078
0.88246 0.80641 0.73843
0.88153 0.80471 0.73608
0.8806 0.80302 0.73374
0.87967 0.80132 0.7314
0.87874 0.79963 0.72906
0.87781 0.79794 0.72672
0.87688 0.79625 0.72439
0.87595 0.79457 0.72206
0.87503 0.79288 0.71974
0.8741 0.7912 0.71742
0.87317 0.78951 0.7151
0.87224 0.78783 0.71279
0.87131 0.78615 0.71047
0.87038 0.78447 0.70817
0.86945 0.78279 0.70586
0.86852 0.78112 0.70356
0.8676 0.77944 0.70126
0.86667 0.77777 0.69897
0.86574 0.7761 0.69667
0.86481 0.77443 0.69439
0.86388 0.77276 0.6921
0.86295 0.77109 0.68982
0.86202 0.76943 0.68754
0.86109 0.76776 0.68527
0.86017 0.7661 0.68299
0.85924 0.76444 0.68073
0.85831 0.76278 0.67846
0.85738 0.76112 0.6762
0.85645 0.75946 0.67394
0.85552 0.7578 0.67169
0.85459 0.75615 0.66943
0.85366 0.75449 0.66719
0.85273 0.75284 0.66494
0.85181 0.75119 0.6627
0.85088 0.74954 0.66046
0.84995 0.7479 0.65822
0.84902 0.74625 0.65599
0.84809 0.74461 0.65376
0.84716 0.74296 0.65154
0.84623 0.74132 0.64932
0.8453 0.73968 0.6471
0.84438 0.73804 0.64488
0.84345 0.73641 0.64267
0.84252 0.73477 0.64046
0.84159 0.73314 0.63826
0.84066 0.7315 0.63606
0.83973 0.72987 0.63386
0.8388 0.72824 0.63166
0.83787 0.72661 0.62947
0.83695 0.72499 0.62728
0.83602 0.72336 0.62509
0.83509 0.72174 0.62291
0.83416 0.72011 0.62073
0.83323 0.71849 0.61856
0.8323 0.71687 0.61639
0.83137 0.71525 0.61422
0.83044 0.71364 0.61205
0.82951 0.71202 0.60989
0.82859 0.71041 0.60773
0.82766 0.7088 0.60558
0.82673 0.70719 0.60342
0.8258 0.70558 0.60127
0.82487 0.70397 0.59913
0.82394 0.70236 0.59699
0.82301 0.70076 0.59485
0.82208 0.69915 0.59271
0.82116 0.69755 0.59058
0.82023 0.69595 0.58845
0.8193 0.69435 0.58632
0.81837 0.69276 0.5842
0.81744 0.69116 0.58208
0.81651 0.68957 0.57997
0.81558 0.68797 0.57785
0.81465 0.68638 0.57575
0.81373 0.68479 0.57364
0.8128 0.6832 0.57154
0.81187 0.68162 0.56944
0.81094 0.68003 0.56734
0.81001 0.67845 0.56525
0.80908 0.67686 0.56316
0.80815 0.67528 0.56107
0.80722 0.6737 0.55899
0.8063 0.67212 0.55691
0.80537 0.67055 0.55484
0.80444 0.66897 0.55276
0.80351 0.6674 0.5507
0.80258 0.66583 0.54863
0.80165 0.66426 0.54657
0.80072 0.66269 0.54451
0.79979 0.66112 0.54245
0.79886 0.65956 0.5404
0.79794 0.65799 0.53835
0.79701 0.65643 0.5363
0.79608 0.65487 0.53426
0.79515 0.65331 0.53222
0.79422 0.65175 0.53018
0.79329 0.65019 0.52815
0.79236 0.64864 0.52612
0.79143 0.64708 0.5241
0.79051 0.64553 0.52207
0.78958 0.64398 0.52005
0.78865 0.64243 0.51804
0.78772 0.64088 0.51602
0.78679 0.63934 0.51402
0.78586 0.63779 0.51201
0.78493 0.63625 0.51001
0.784 0.63471 0.50801
0.78308 0.63317 0.50601
0.78215 0.63163 0.50402
0.78122 0.63009 0.50203
0.78029 0.62856 0.50004
0.77936 0.62702 0.49806
0.77843 0.62549 0.49608
0.77743 0.6237 0.49365
0.77644 0.62192 0.49123
0.77544 0.62014 0.48881
0.77444 0.61836 0.48639
0.77344 0.61658 0.48399
0.77245 0.61481 0.48158
0.77145 0.61304 0.47918
0.77045 0.61127 0.47679
0.76945 0.6095 0.47439
0.76846 0.60774 0.47201
0.76746 0.60597 0.46963
0.76646 0.60421 0.46725
0.76546 0.60246 0.46488
0.76447 0.6007 0.46251
0.76347 0.59895 0.46015
0.76247 0.5972 0.45779
0.76147 0.59545 0.45543
0.76047 0.5937 0.45308
0.75948 0.59196 0.45074
0.75848 0.59022 0.4484
0.75748 0.58848 0.44606
0.75648 0.58674 0.44373
0.75549 0.58501 0.44141
0.75449 0.58328 0.43908
0.75349 0.58155 0.43677
0.75249 0.57982 0.43445
0.7515 0.5781 0.43215
0.7505 0.57638 0.42984
0.7495 0.57466 0.42754
0.7485 0.57294 0.42525
0.74751 0.57123 0.42296
0.74651 0.56951 0.42068
0.74551 0.5678 0.4184
0.74451 0.56609 0.41612
0.74352 0.56439 0.41385
0.74252 0.56269 0.41158
0.74152 0.56099 0.40932
0.74052 0.55929 0.40706
0.73953 0.55759 0.40481
0.73853 0.5559 0.40256
0.73753 0.55421 0.40032
0.73653 0.55252 0.39808
0.73553 0.55083 0.39585
0.73454 0.54915 0.39362
0.73354 0.54747 0.39139
0.73254 0.54579 0.38917
0.73154 0.54411 0.38696
0.73055 0.54244 0.38475
0.72955 0.54077 0.38254
0.72855 0.5391 0.38034
0.72755 0.53743 0.37814
0.72656 0.53576 0.37595
0.72556 0.5341 0.37376
0.72456 0.53244 0.37158
0.72356 0.53078 0.3694
0.72257 0.52913 0.36722
0.72157 0.52748 0.36505
0.72057 0.52583 0.36289
0.71957 0.52418 0.36073
0.71858 0.52253 0.35857
0.71758 0.52089 0.35642
0.71658 0.51925 0.35428
0.71558 0.51761 0.35213
0.71459 0.51597 0.35
0.71359 0.51434 0.34787
0.71259 0.51271 0.34574
0.71159 0.51108 0.34361
0.7106 0.50945 0.3415
0.7096 0.50783 0.33938
0.7086 0.50621 0.33727
0.7076 0.50459 0.33517
0.7066 0.50297 0.33307
0.70561 0.50136 0.33097
0.70461 0.49975 0.32888
0.70361 0.49814 0.32679
0.70261 0.49653 0.32471
0.70162 0.49493 0.32263
0.70062 0.49332 0.32056
0.69962 0.49172 0.31849
0.69862 0.49013 0.31643
0.69763 0.48853 0.31437
0.69663 0.48694 0.31232
0.69563 0.48535 0.31027
0.69463 0.48376 0.30822
0.69364 0.48218 0.30618
0.69264 0.48059 0.30415
0.69164 0.47901 0.30212
0.69064 0.47743 0.30009
0.68965 0.47586 0.29807
0.68865 0.47429 0.29605
0.68765 0.47271 0.29404
0.68665 0.47115 0.29203
0.68566 0.46958 0.29003
0.68466 0.46802 0.28803
0.68366 0.46646 0.28603
0.68266 0.4649 0.28404
0.68166 0.46334 0.28206
0.68067 0.46179 0.28008
0.67967 0.46023 0.2781
0.67867 0.45869 0.27613
0.67767 0.45714 0.27416
0.67668 0.4556 0.2722
0.67568 0.45405 0.27025
0.67468 0.45251 0.26829
0.67368 0.45098 0.26634
0.67269 0.44944 0.2644
0.67169 0.44791 0.26246
0.67069 0.44638 0.26053
0.66969 0.44485 0.2586
0.6687 0.44333 0.25667
0.6677 0.44181 0.25475
0.6667 0.44029 0.25284
0.6657 0.43877 0.25093
0.66471 0.43725 0.24902
0.66373 0.43538 0.24652
0.66276 0.43351 0.24402
0.66179 0.43164 0.24153
0.66082 0.42978 0.23904
0.65985 0.42792 0.23656
0.65888 0.42606 0.23409
0.6579 0.4242 0.23162
0.65693 0.42235 0.22916
0.65596 0.4205 0.22671
0.65499 0.41866 0.22426
0.65402 0.41681 0.22182
0.65305 0.41497 0.21939
0.65207 0.41314 0.21696
0.6511 0.4113 0.21454
0.65013 0.40947 0.21212
0.64916 0.40764 0.20971
0.64819 0.40582 0.20731
0.64722 0.40399 0.20491
0.64624 0.40218 0.20252
0.64527 0.40036 0.20014
0.6443 0.39855 0.19776
0.64333 0.39674 0.19539
0.64236 0.39493 0.19302
0.64138 0.39312 0.19066
0.64041 0.39132 0.18831
0.63944 0.38953 0.18596
0.63847 0.38773 0.18362
0.6375 0.38594 0.18129
0.63653 0.38415 0.17896
0.63555 0.38236 0.17664
0.63458 0.38058 0.17432
0.63361 0.3788 0.17202
0.63264 0.37703 0.16971
0.63167 0.37525 0.16742
0.6307 0.37348 0.16513
0.62972 0.37172 0.16284
0.62875 0.36995 0.16057
0.62778 0.36819 0.15829
0.62681 0.36643 0.15603
0.62584 0.36468 0.15377
0.62487 0.36292 0.15152
0.62389 0.36118 0.14927
0.62292 0.35943 0.14703
0.62195 0.35769 0.1448
0.62098 0.35595 0.14257
0.62001 0.35421 0.14035
0.61904 0.35248 0.13814
0.61806 0.35075 0.13593
0.61709 0.34902 0.13373
0.61612 0.3473 0.13153
0.61515 0.34558 0.12934
0.61418 0.34386 0.12716
0.6132 0.34214 0.12498
0.61223 0.34043 0.12281
0.61126 0.33872 0.12065
0.61029 0.33702 0.11849
0.60932 0.33532 0.11634
0.60835 0.33362 0.11419
0.60737 0.33192 0.11205
0.6064 0.33023 0.10992
0.60543 0.32854 0.10779
0.60446 0.32685 0.10567
0.60349 0.32517 0.10356
0.60252 0.32349 0.10145
0.60154 0.32181 0.099351
0.60057 0.32014 0.097256
0.5996 0.31847 0.095167
0.59863 0.3168 0.093084
0.59766 0.31514 0.091008
0.59669 0.31348 0.088938
0.59571 0.31182 0.086875
0.59474 0.31017 0.084817
0.59377 0.30851 0.082766
0.5928 0.30687 0.080722
0.59183 0.30522 0.078683
0.59086 0.30358 0.076651
0.58988 0.30194 0.074625
0.58891 0.3003 0.072605
0.58794 0.29867 0.070592
0.58697 0.29704 0.068585
0.586 0.29542 0.066584
0.58503 0.29379 0.06459
0.58405 0.29218 0.062602
0.58308 0.29056 0.06062
0.58211 0.28895 0.058644
0.58114 0.28734 0.056675
0.58017 0.28573 0.054712
0.57919 0.28413 0.052756
0.57822 0.28253 0.050805
0.57725 0.28093 0.048861
0.57628 0.27934 0.046923
0.57531 0.27774 0.044992
0.57434 0.27616 0.043066
0.57336 0.27457 0.041147
0.57239 0.27299 0.039235
0.57142 0.27142 0.037328
0.57045 0.26984 0.035428
0.56948 0.26827 0.033535
0.56851 0.2667 0.031647
0.56753 0.26514 0.029766
0.56656 0.26358 0.027891
0.56559 0.26202 0.026022
0.56462 0.26046 0.02416
0.56365 0.25891 0.022304
0.56268 0.25736 0.020454
0.5617 0.25582 0.018611
0.56073 0.25428 0.016773
0.55976 0.25274 0.014943
0.55879 0.25121 0.013118
0.55782 0.24967 0.0113
0.55685 0.24814 0.0094876
0.55587 0.24662 0.0076818
0.5549 0.2451 0.0058824
0.55394 0.2428 0.005884
0.55298 0.24051 0.0058857
0.55201 0.23822 0.0058873
0.55105 0.23594 0.0058889
0.55009 0.23367 0.0058904
0.54912 0.2314 0.005892
0.54816 0.22915 0.0058935
0.5472 0.22689 0.0058949
0.54623 0.22465 0.0058964
0.54527 0.22241 0.0058978
0.54431 0.22017 0.0058992
0.54334 0.21795 0.0059005
0.54238 0.21573 0.0059018
0.54142 0.21352 0.0059031
0.54045 0.21131 0.0059044
0.53949 0.20911 0.0059056
0.53853 0.20692 0.0059068
0.53756 0.20473 0.005908
0.5366 0.20255 0.0059092
0.53564 0.20038 0.0059103
0.53467 0.19821 0.0059114
0.53371 0.19605 0.0059124
0.53275 0.1939 0.0059134
0.53179 0.19175 0.0059144
0.53082 0.18961 0.0059154
0.52986 0.18747 0.0059163
0.5289 0.18535 0.0059172
0.52793 0.18323 0.0059181
0.52697 0.18111 0.005919
0.52601 0.179 0.0059198
0.52504 0.1769 0.0059206
0.52408 0.17481 0.0059213
0.52312 0.17272 0.005922
0.52215 0.17064 0.0059227
0.52119 0.16857 0.0059234
0.52023 0.1665 0.005924
0.51926 0.16444 0.0059247
0.5183 0.16238 0.0059252
0.51734 0.16033 0.0059258
0.51637 0.15829 0.0059263
0.51541 0.15626 0.0059268
0.51445 0.15423 0.0059273
0.51348 0.15221 0.0059277
0.51252 0.15019 0.0059281
0.51156 0.14818 0.0059285
0.5106 0.14618 0.0059288
0.50963 0.14419 0.0059291
0.50867 0.1422 0.0059294
0.50771 0.14022 0.0059296
0.50674 0.13824 0.0059299
0.50578 0.13627 0.0059301
0.50482 0.13431 0.0059302
0.50385 0.13235 0.0059304
0.50289 0.1304 0.0059305
0.50193 0.12846 0.0059305
0.50096 0.12652 0.0059306
0.5 0.1246 0.0059306
0.49904 0.12267 0.0059306
0.49807 0.12076 0.0059305
0.49711 0.11885 0.0059305
0.49615 0.11694 0.0059304
0.49518 0.11505 0.0059302
0.49422 0.11316 0.0059301
0.49326 0.11127 0.0059299
0.49229 0.1094 0.0059296
0.49133 0.10753 0.0059294
0.49037 0.10566 0.0059291
0.4894 0.10381 0.0059288
0.48844 0.10195 0.0059285
0.48748 0.10011 0.0059281
0.48652 0.098274 0.0059277
0.48555 0.096443 0.0059273
0.48459 0.094619 0.0059268
0.48363 0.092801 0.0059263
0.48266 0.09099 0.0059258
0.4817 0.089186 0.0059252
0.48074 0.087388 0.0059247
0.47977 0.085597 0.005924
0.47881 0.083812 0.0059234
0.47785 0.082034 0.0059227
0.47688 0.080262 0.005922
0.47592 0.078497 0.0059213
0.47496 0.076739 0.0059206
0.47399 0.074987 0.0059198
0.47303 0.073242 0.005919
0.47207 0.071503 0.0059181
0.4711 0.069771 0.0059172
0.47014 0.068046 0.0059163
0.46918 0.066327 0.0059154
0.46821 0.064614 0.0059144
0.46725 0.062909 0.0059134
0.46629 0.061209 0.0059124
0.46533 0.059517 0.0059114
0.46436 0.057831 0.0059103
0.4634 0.056151 0.0059092
0.46244 0.054478 0.005908
0.46147 0.052812 0.0059068
0.46051 0.051152 0.0059056
0.45955 0.049499 0.0059044
0.45858 0.047853 0.0059031
0.45762 0.046213 0.0059018
0.45666 0.044579 0.0059005
0.45569 0.042952 0.0058992
0.45473 0.041332 0.0058978
0.45377 0.039718 0.0058964
0.4528 0.038111 0.0058949
0.45184 0.03651 0.0058935
0.45088 0.034916 0.005892
0.44991 0.033329 0.0058904
0.44895 0.031748 0.0058889
0.44799 0.030174 0.0058873
0.44702 0.028606 0.0058857
0.44606 0.027045 0.005884
0.4451 0.02549 0.0058824
0.44413 0.025235 0.0058395
0.44317 0.02498 0.0057968
0.44221 0.024726 0.0057543
0.44125 0.024473 0.0057119
0.44028 0.024221 0.0056696
0.43932 0.02397 0.0056275
0.43836 0.023719 0.0055856
0.43739 0.02347 0.0055438
0.43643 0.023221 0.0055021
0.43547 0.022973 0.0054606
0.4345 0.022726 0.0054192
0.43354 0.02248 0.005378
0.43258 0.022235 0.005337
0.43161 0.021991 0.0052961
0.43065 0.021747 0.0052553
0.42969 0.021505 0.0052147
0.42872 0.021263 0.0051743
0.42776 0.021022 0.005134
0.4268 0.020782 0.0050938
0.42583 0.020543 0.0050538
0.42487 0.020305 0.0050139
0.42391 0.020067 0.0049742
0.42294 0.019831 0.0049347
0.42198 0.019595 0.0048953
0.42102 0.01936 0.004856
0.42006 0.019126 0.0048169
0.41909 0.018893 0.0047779
0.41813 0.018661 0.0047391
0.41717 0.01843 0.0047005
0.4162 0.018199 0.0046619
0.41524 0.01797 0.0046236
0.41428 0.017741 0.0045854
0.41331 0.017513 0.0045473
0.41235 0.017286 0.0045094
0.41139 0.01706 0.0044716
0.41042 0.016835 0.004434
0.40946 0.016611 0.0043965
0.4085 0.016387 0.0043592
0.40753 0.016165 0.0043221
0.40657 0.015943 0.004285
0.40561 0.015722 0.0042482
0.40464 0.015502 0.0042115
0.40368 0.015283 0.0041749
0.40272 0.015064 0.0041385
0.40175 0.014847 0.0041022
0.40079 0.014631 0.0040661
0.39983 0.014415 0.0040301
0.39886 0.0142 0.0039943
0.3979 0.013986 0.0039586
0.39694 0.013773 0.0039231
0.39598 0.013561 0.0038877
0.39501 0.01335 0.0038525
0.39405 0.013139 0.0038174
0.39309 0.01293 0.0037825
0.39212 0.012721 0.0037477
0.39116 0.012513 0.0037131
0.3902 0.012306 0.0036786
0.38923 0.0121 0.0036443
0.38827 0.011895 0.0036101
0.38731 0.011691 0.0035761
0.38634 0.011487 0.0035422
0.38538 0.011285 0.0035085
0.38442 0.011083 0.0034749
0.38345 0.010882 0.0034415
0.38249 0.010682 0.0034082
0.38153 0.010483 0.0033751
0.38056 0.010285 0.0033421
0.3796 0.010088 0.0033093
0.37864 0.0098913 0.0032766
0.37767 0.0096957 0.0032441
0.37671 0.009501 0.0032117
0.37575 0.0093072 0.0031795
0.37479 0.0091143 0.0031474
0.37382 0.0089223 0.0031155
0.37286 0.0087311 0.0030837
0.3719 0.0085408 0.003052
0.37093 0.0083514 0.0030206
0.36997 0.0081628 0.0029892
0.36901 0.0079751 0.002958
0.36804 0.0077883 0.002927
0.36708 0.0076024 0.0028961
0.36612 0.0074174 0.0028654
0.36515 0.0072332 0.0028348
0.36419 0.00705 0.0028044
0.36323 0.0068676 0.0027741
0.36226 0.006686 0.0027439
0.3613 0.0065054 0.0027139
0.36034 0.0063256 0.0026841
0.35937 0.0061467 0.0026544
0.35841 0.0059687 0.0026249
0.35745 0.0057916 0.0025955
0.35648 0.0056153 0.0025662
0.35552 0.00544 0.0025372
0.35456 0.0052655 0.0025082
0.35359 0.0050919 0.0024794
0.35263 0.0049192 0.0024508
0.35167 0.0047473 0.0024223
0.35071 0.0045763 0.0023939
0.34974 0.0044063 0.0023657
0.34878 0.0042371 0.0023377
0.34782 0.0040687 0.0023098
0.34685 0.0039013 0.0022821
0.34589 0.0037347 0.0022545
0.34493 0.003569 0.002227
0.34396 0.0034042 0.0021997
0.343 0.0032403 0.0021726
0.34204 0.0030773 0.0021456
0.34107 0.0029151 0.0021187
0.34011 0.0027539 0.002092
0.33915 0.0025935 0.0020655
0.33818 0.002434 0.0020391
0.33722 0.0022754 0.0020128
0.33626 0.0021176 0.0019867
0.33529 0.0019608 0.0019608
0.33432 0.0019608 0.0019608
0.33335 0.0019608 0.0019608
0.33238 0.0019608 0.0019608
0.33141 0.0019608 0.0019608
0.33044 0.0019608 0.0019608
0.32946 0.0019608 0.0019608
0.32849 0.0019608 0.0019608
0.32752 0.0019608 0.0019608
0.32655 0.0019608 0.0019608
0.32558 0.0019608 0.0019608
0.32461 0.0019608 0.0019608
0.32363 0.0019608 0.0019608
0.32266 0.0019608 0.0019608
0.32169 0.0019608 0.0019608
0.32072 0.0019608 0.0019608
0.31975 0.0019608 0.0019608
0.31877 0.0019608 0.0019608
0.3178 0.0019608 0.0019608
0.31683 0.0019608 0.0019608
0.31586 0.0019608 0.0019608
0.31489 0.0019608 0.0019608
0.31392 0.0019608 0.0019608
0.31294 0.0019608 0.0019608
0.31197 0.0019608 0.0019608
0.311 0.0019608 0.0019608
0.31003 0.0019608 0.0019608
0.30906 0.0019608 0.0019608
0.30809 0.0019608 0.0019608
0.30711 0.0019608 0.0019608
0.30614 0.0019608 0.0019608
0.30517 0.0019608 0.0019608
0.3042 0.0019608 0.0019608
0.30323 0.0019608 0.0019608
0.30226 0.0019608 0.0019608
0.30128 0.0019608 0.0019608
0.30031 0.0019608 0.0019608
0.29934 0.0019608 0.0019608
0.29837 0.0019608 0.0019608
0.2974 0.0019608 0.0019608
0.29643 0.0019608 0.0019608
0.29545 0.0019608 0.0019608
0.29448 0.0019608 0.0019608
0.29351 0.0019608 0.0019608
0.29254 0.0019608 0.0019608
0.29157 0.0019608 0.0019608
0.2906 0.0019608 0.0019608
0.28962 0.0019608 0.0019608
0.28865 0.0019608 0.0019608
0.28768 0.0019608 0.0019608
0.28671 0.0019608 0.0019608
0.28574 0.0019608 0.0019608
0.28476 0.0019608 0.0019608
0.28379 0.0019608 0.0019608
0.28282 0.0019608 0.0019608
0.28185 0.0019608 0.0019608
0.28088 0.0019608 0.0019608
0.27991 0.0019608 0.0019608
0.27893 0.0019608 0.0019608
0.27796 0.0019608 0.0019608
0.27699 0.0019608 0.0019608
0.27602 0.0019608 0.0019608
0.27505 0.0019608 0.0019608
0.27408 0.0019608 0.0019608
0.2731 0.0019608 0.0019608
0.27213 0.0019608 0.0019608
0.27116 0.0019608 0.0019608
0.27019 0.0019608 0.0019608
0.26922 0.0019608 0.0019608
0.26825 0.0019608 0.0019608
0.26727 0.0019608 0.0019608
0.2663 0.0019608 0.0019608
0.26533 0.0019608 0.0019608
0.26436 0.0019608 0.0019608
0.26339 0.0019608 0.0019608
0.26242 0.0019608 0.0019608
0.26144 0.0019608 0.0019608
0.26047 0.0019608 0.0019608
0.2595 0.0019608 0.0019608
0.25853 0.0019608 0.0019608
0.25756 0.0019608 0.0019608
0.25659 0.0019608 0.0019608
0.25561 0.0019608 0.0019608
0.25464 0.0019608 0.0019608
0.25367 0.0019608 0.0019608
0.2527 0.0019608 0.0019608
0.25173 0.0019608 0.0019608
0.25075 0.0019608 0.0019608
0.24978 0.0019608 0.0019608
0.24881 0.0019608 0.0019608
0.24784 0.0019608 0.0019608
0.24687 0.0019608 0.0019608
0.2459 0.0019608 0.0019608
0.24492 0.0019608 0.0019608
0.24395 0.0019608 0.0019608
0.24298 0.0019608 0.0019608
0.24201 0.0019608 0.0019608
0.24104 0.0019608 0.0019608
0.24007 0.0019608 0.0019608
0.23909 0.0019608 0.0019608
0.23812 0.0019608 0.0019608
0.23715 0.0019608 0.0019608
0.23618 0.0019608 0.0019608
0.23521 0.0019608 0.0019608
0.23424 0.0019608 0.0019608
0.23326 0.0019608 0.0019608
0.23229 0.0019608 0.0019608
0.23132 0.0019608 0.0019608
0.23035 0.0019608 0.0019608
0.22938 0.0019608 0.0019608
0.22841 0.0019608 0.0019608
0.22743 0.0019608 0.0019608
0.22646 0.0019608 0.0019608
0.22549 0.0019608 0.0019608
0.22456 0.0019608 0.0019608
0.22363 0.0019608 0.0019608
0.2227 0.0019608 0.0019608
0.22178 0.0019608 0.0019608
0.22085 0.0019608 0.0019608
0.21992 0.0019608 0.0019608
0.21899 0.0019608 0.0019608
0.21806 0.0019608 0.0019608
0.21713 0.0019608 0.0019608
0.2162 0.0019608 0.0019608
0.21527 0.0019608 0.0019608
0.21434 0.0019608 0.0019608
0.21342 0.0019608 0.0019608
0.21249 0.0019608 0.0019608
0.21156 0.0019608 0.0019608
0.21063 0.0019608 0.0019608
0.2097 0.0019608 0.0019608
0.20877 0.0019608 0.0019608
0.20784 0.0019608 0.0019608
0.20691 0.0019608 0.0019608
0.20599 0.0019608 0.0019608
0.20506 0.0019608 0.0019608
0.20413 0.0019608 0.0019608
0.2032 0.0019608 0.0019608
0.20227 0.0019608 0.0019608
0.20134 0.0019608 0.0019608
0.20041 0.0019608 0.0019608
0.19948 0.0019608 0.0019608
0.19856 0.0019608 0.0019608
0.19763 0.0019608 0.0019608
0.1967 0.0019608 0.0019608
0.19577 0.0019608 0.0019608
0.19484 0.0019608 0.0019608
0.19391 0.0019608 0.0019608
0.19298 0.0019608 0.0019608
0.19205 0.0019608 0.0019608
0.19112 0.0019608 0.0019608
0.1902 0.0019608 0.0019608
0.18927 0.0019608 0.0019608
0.18834 0.0019608 0.0019608
0.18741 0.0019608 0.0019608
0.18648 0.0019608 0.0019608
0.18555 0.0019608 0.0019608
0.18462 0.0019608 0.0019608
0.18369 0.0019608 0.0019608
0.18277 0.0019608 0.0019608
0.18184 0.0019608 0.0019608
0.18091 0.0019608 0.0019608
0.17998 0.0019608 0.0019608
0.17905 0.0019608 0.0019608
0.17812 0.0019608 0.0019608
0.17719 0.0019608 0.0019608
0.17626 0.0019608 0.0019608
0.17534 0.0019608 0.0019608
0.17441 0.0019608 0.0019608
0.17348 0.0019608 0.0019608
0.17255 0.0019608 0.0019608
0.17162 0.0019608 0.0019608
0.17069 0.0019608 0.0019608
0.16976 0.0019608 0.0019608
0.16883 0.0019608 0.0019608
0.16791 0.0019608 0.0019608
0.16698 0.0019608 0.0019608
0.16605 0.0019608 0.0019608
0.16512 0.0019608 0.0019608
0.16419 0.0019608 0.0019608
0.16326 0.0019608 0.0019608
0.16233 0.0019608 0.0019608
0.1614 0.0019608 0.0019608
0.16047 0.0019608 0.0019608
0.15955 0.0019608 0.0019608
0.15862 0.0019608 0.0019608
0.15769 0.0019608 0.0019608
0.15676 0.0019608 0.0019608
0.15583 0.0019608 0.0019608
0.1549 0.0019608 0.0019608
0.15397 0.0019608 0.0019608
0.15304 0.0019608 0.0019608
0.15212 0.0019608 0.0019608
0.15119 0.0019608 0.0019608
0.15026 0.0019608 0.0019608
0.14933 0.0019608 0.0019608
0.1484 0.0019608 0.0019608
0.14747 0.0019608 0.0019608
0.14654 0.0019608 0.0019608
0.14561 0.0019608 0.0019608
0.14469 0.0019608 0.0019608
0.14376 0.0019608 0.0019608
0.14283 0.0019608 0.0019608
0.1419 0.0019608 0.0019608
0.14097 0.0019608 0.0019608
0.14004 0.0019608 0.0019608
0.13911 0.0019608 0.0019608
0.13818 0.0019608 0.0019608
0.13725 0.0019608 0.0019608
0.13633 0.0019608 0.0019608
0.1354 0.0019608 0.0019608
0.13447 0.0019608 0.0019608
0.13354 0.0019608 0.0019608
0.13261 0.0019608 0.0019608
0.13168 0.0019608 0.0019608
0.13075 0.0019608 0.0019608
0.12982 0.0019608 0.0019608
0.1289 0.0019608 0.0019608
0.12797 0.0019608 0.0019608
0.12704 0.0019608 0.0019608
0.12611 0.0019608 0.0019608
0.12518 0.0019608 0.0019608
0.12425 0.0019608 0.0019608
0.12332 0.0019608 0.0019608
0.12239 0.0019608 0.0019608
0.12147 0.0019608 0.0019608
0.12054 0.0019608 0.0019608
0.11961 0.0019608 0.0019608
0.11858 0.0029837 0.0029837
0.11754 0.0039886 0.0039886
0.11651 0.0049753 0.0049753
0.11548 0.0059439 0.0059439
0.11445 0.0068944 0.0068944
0.11342 0.0078268 0.0078268
0.11238 0.0087412 0.0087412
0.11135 0.0096374 0.0096374
0.11032 0.010515 0.010515
0.10929 0.011375 0.011375
0.10826 0.012217 0.012217
0.10722 0.013041 0.013041
0.10619 0.013847 0.013847
0.10516 0.014634 0.014634
0.10413 0.015404 0.015404
0.1031 0.016155 0.016155
0.10206 0.016888 0.016888
0.10103 0.017604 0.017604
0.1 0.018301 0.018301
0.098968 0.01898 0.01898
0.097936 0.01964 0.01964
0.096904 0.020283 0.020283
0.095872 0.020908 0.020908
0.09484 0.021514 0.021514
0.093808 0.022103 0.022103
0.092776 0.022673 0.022673
0.091744 0.023225 0.023225
0.090712 0.023759 0.023759
0.08968 0.024275 0.024275
0.088648 0.024773 0.024773
0.087616 0.025253 0.025253
0.086584 0.025715 0.025715
0.085552 0.026158 0.026158
0.08452 0.026584 0.026584
0.083488 0.026991 0.026991
0.082456 0.02738 0.02738
0.081424 0.027752 0.027752
0.080392 0.028105 0.028105
0.07936 0.02844 0.02844
0.078328 0.028756 0.028756
0.077296 0.029055 0.029055
0.076264 0.029336 0.029336
0.075232 0.029598 0.029598
0.0742 0.029843 0.029843
0.073168 0.030069 0.030069
0.072136 0.030277 0.030277
0.071104 0.030467 0.030467
0.070072 0.030639 0.030639
0.06904 0.030793 0.030793
0.068008 0.030929 0.030929
0.066976 0.031047 0.031047
0.065944 0.031146 0.031146
0.064912 0.031228 0.031228
0.06388 0.031291 0.031291
0.062848 0.031336 0.031336
0.061816 0.031363 0.031363
0.060784 0.031373 0.031373
0.059752 0.031363 0.031363
0.05872 0.031336 0.031336
0.057688 0.031291 0.031291
0.056656 0.031228 0.031228
0.055624 0.031146 0.031146
0.054592 0.031047 0.031047
0.05356 0.030929 0.030929
0.052528 0.030793 0.030793
0.051496 0.030639 0.030639
0.050464 0.030467 0.030467
0.049432 0.030277 0.030277
0.0484 0.030069 0.030069
0.047368 0.029843 0.029843
0.046336 0.029598 0.029598
0.045304 0.029336 0.029336
0.044272 0.029055 0.029055
0.04324 0.028756 0.028756
0.042208 0.02844 0.02844
0.041176 0.028105 0.028105
0.040144 0.027752 0.027752
0.039112 0.02738 0.02738
0.03808 0.026991 0.026991
0.037049 0.026584 0.026584
0.036017 0.026158 0.026158
0.034985 0.025715 0.025715
0.033953 0.025253 0.025253
0.032921 0.024773 0.024773
0.031889 0.024275 0.024275
0.030857 0.023759 0.023759
0.029825 0.023225 0.023225
0.028793 0.022673 0.022673
0.027761 0.022103 0.022103
0.026729 0.021514 0.021514
0.025697 0.020908 0.020908
0.024665 0.020283 0.020283
0.023633 0.01964 0.01964
0.022601 0.01898 0.01898
0.021569 0.018301 0.018301
0.020537 0.017604 0.017604
0.019505 0.016888 0.016888
0.018473 0.016155 0.016155
0.017441 0.015404 0.015404
0.016409 0.014634 0.014634
0.015377 0.013847 0.013847
0.014345 0.013041 0.013041
0.013313 0.012217 0.012217
0.012281 0.011375 0.011375
0.011249 0.010515 0.010515
0.010217 0.0096374 0.0096374
0.0091847 0.0087412 0.0087412
0.0081527 0.0078268 0.0078268
0.0071207 0.0068944 0.0068944
0.0060888 0.0059439 0.0059439
0.0050568 0.0049753 0.0049753
0.0040248 0.0039886 0.0039886
0.0029928 0.0029837 0.0029837
0.0019608 0.0019608 0.0019608
];


if nargin==1
    a=interp1(1:1024,phasemap(:,1),1:1024/npx:1024).';
    b=interp1(1:1024,phasemap(:,2),1:1024/npx:1024).';
    c=interp1(1:1024,phasemap(:,3),1:1024/npx:1024).';
    phasemap=[a b c];
end




























