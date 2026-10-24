# Saved eight-atom counts, independently parsed

| Operator | Target <- source | Offsets | Shell range | Related six-atom block/count |
|---|---|---:|---|---|
| sigma2P | Om <- Om | 571 | [-7, 7] | (U,U) / 571 |
| sigma2P | P <- P | 1195 | [-9, 9] | (V,V) / 1617 |
| sigma2P | P <- lb | 709 | [-8, 8] | (V,U) / 1357 |
| sigma2P | al <- Om | 461 | [-6, 6] | (x,U) / 461 |
| sigma2P | al <- al | 571 | [-7, 7] | (x,x) / 571 |
| sigma2P | beP <- Om | 709 | [-8, 8] | (V,U) / 1357 |
| sigma2P | beP <- beP | 1195 | [-9, 9] | (V,V) / 1617 |
| sigma2P | lR <- lR | 571 | [-7, 7] | (Slog,Slog) / 571 |
| sigma2P | la <- la | 571 | [-7, 7] | (L,L) / 571 |
| sigma2P | lb <- lb | 571 | [-7, 7] | (U,U) / 571 |
| sigma2P | one <- al | 971 | [-9, 9] | (one,x) / 1321 |
| sigma2P | one <- beP | 1371 | [-10, 10] | (one,V) / 1371 |
| sigma2P | one <- lR | 915 | [-9, 9] | (one,Slog) / 1237 |
| sigma2P | one <- la | 939 | [-9, 9] | (one,L) / 1299 |
| sigma2P | one <- lb | 1037 | [-9, 9] | (one,U) / 1407 |
| sigma2P | one <- one | 1209 | [-9, 9] | (one,one) / 1617 |
| C1 | Om <- Om | 419 | [-7, 7] | (U,U) / 419 |
| C1 | P <- P | 947 | [-9, 9] | (V,V) / 1297 |
| C1 | P <- lb | 369 | [-7, 7] | (V,U) / 897 |
| C1 | al <- al | 419 | [-7, 7] | (x,x) / 419 |
| C1 | beP <- Om | 369 | [-7, 7] | (V,U) / 897 |
| C1 | beP <- beP | 947 | [-9, 9] | (V,V) / 1297 |
| C1 | lR <- lR | 419 | [-7, 7] | (Slog,Slog) / 419 |
| C1 | la <- la | 419 | [-7, 7] | (L,L) / 419 |
| C1 | lb <- lb | 419 | [-7, 7] | (U,U) / 419 |
| C1 | one <- al | 625 | [-8, 8] | (one,x) / 861 |
| C1 | one <- beP | 931 | [-9, 9] | (one,V) / 931 |
| C1 | one <- lR | 551 | [-8, 8] | (one,Slog) / 727 |
| C1 | one <- la | 617 | [-8, 8] | (one,L) / 817 |
| C1 | one <- lb | 667 | [-8, 8] | (one,U) / 963 |
| C1 | one <- one | 967 | [-9, 9] | (one,one) / 1309 |
| C2 | Om <- Om | 251 | [-6, 6] | (U,U) / 251 |
| C2 | P <- P | 683 | [-8, 8] | (V,V) / 931 |
| C2 | al <- al | 251 | [-6, 6] | (x,x) / 251 |
| C2 | beP <- beP | 683 | [-8, 8] | (V,V) / 931 |
| C2 | lR <- lR | 251 | [-6, 6] | (Slog,Slog) / 251 |
| C2 | la <- la | 251 | [-6, 6] | (L,L) / 251 |
| C2 | lb <- lb | 251 | [-6, 6] | (U,U) / 251 |
| C2 | one <- one | 683 | [-8, 8] | (one,one) / 955 |

Related means only basis correspondence; U=Om+pi lb and V=(beP+pi P/2)/r12 mix blocks and derivatives. Compare after removing row factors.

Row factors: `{'Om': 'u*v*w', 'al': 'u*v*w', 'beP': 'u*v*w*(u**2 - 2*u*v + 2*u*w + v**2 + 2*v*w + w**2)', 'lR': 'u*v*w', 'la': 'u*v*w', 'lb': 'u*v*w', 'one': 'u*v*w*(2*u**2 + 2*u*w + 2*v**2 + 2*v*w + w**2)/2', 'P': 'u*v*w*(u**2 - 2*u*v + 2*u*w + v**2 + 2*v*w + w**2)'}`.
