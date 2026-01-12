# SPRAWOZDANIE Z LABORATORIUM NR 5

**Imię i nazwisko:** Jakub Sputo, Jakub Stanula-Kaczka, Kacper Strzesak  
**Kierunek:** Informatyka techniczna  
**Kurs:** Optymalizacja  
**Grupa:** 3  
**Data zajęć:** [Data do uzupełnienia]  
**Numer ćwiczenia:** 5  
**Temat ćwiczenia:** Optymalizacja wielokryterialna  

---

## Cel ćwiczenia

Celem ćwiczenia było zapoznanie się z problematyką optymalizacji wielokryterialnej oraz wyznaczenie rozwiązań minimalnych w sensie Pareto. Ćwiczenie obejmowało implementację i praktyczne zastosowanie metody Powella w połączeniu z metodą kryterium ważonego do rozwiązywania problemów wielokryterialnych.

---

## Opis zadania

Zadanie zostało podzielone na dwie części:

### Część 1: Testowa funkcja celu

Zaimplementowano optymalizację dwukryterialną z funkcjami celu danymi wzorami:

$$f_1(x_1, x_2) = a \cdot ((x_1 - 3)^2 + (x_2 - 3)^2)$$

$$f_2(x_1, x_2) = \frac{1}{a} \cdot ((x_1 + 3)^2 + (x_2 + 3)^2)$$

gdzie parametr $a$ przyjmował wartości: 1, 10, 100.

Przeprowadzono 101 optymalizacji dla każdej wartości parametru $a$, przy współczynniku wagowym $w \in \{0, 0.01, 0.02, ..., 1.0\}$. Punkty startowe wybierano losowo z przedziału $x_1^{(0)}, x_2^{(0)} \in [-10, 10]$. 

Dla zapewnienia uczciwego porównania między różnymi wartościami parametru $a$, dla każdego $w$ wszystkie trzy wartości $a$ wykorzystywały identyczny punkt startowy.

### Część 2: Problem rzeczywisty - belka wspornikowa

Rozwiązano problem rzeczywisty dotyczący optymalizacji belki wspornikowej o przekroju kołowym. Belka o długości $l$ i średnicy $d$ jest obciążona siłą $P = 2$ kN.

**Funkcje celu:**
- $f_1$ - masa belki: $m = \rho \cdot \pi \cdot \left(\frac{d}{2}\right)^2 \cdot l$
- $f_2$ - ugięcie belki: $u = \frac{64 \cdot P \cdot l^3}{3 \cdot E \cdot \pi \cdot d^4}$

**Parametry materiałowe:**
- $\rho = 8920$ kg/m³ (gęstość)
- $E = 120$ GPa (moduł Younga)

**Ograniczenia:**
- Ugięcie: $u \leq 2.5$ mm
- Naprężenie: $\sigma = \frac{32 \cdot P \cdot l}{\pi \cdot d^3} \leq 300$ MPa
- Zakres zmiennych: $l \in [200, 1000]$ mm, $d \in [10, 50]$ mm

Przeprowadzono 101 optymalizacji dla współczynnika wagowego $w \in \{0, 0.01, 0.02, ..., 1.0\}$ z losowymi punktami startowymi.

---

## Algorytmy optymalizacji

### Metoda kryterium ważonego

Problem wielokryterialny zamieniono na problem jednokryterialny stosując metodę kryterium ważonego:

$$f(x) = w \cdot f_1(x) + (1 - w) \cdot f_2(x)$$

gdzie $w \in [0, 1]$ jest współczynnikiem wagowym. Zmiana wartości $w$ pozwala na wyznaczenie całego frontu Pareto.

### Metoda Powella

Do wyznaczenia minimum funkcji celu zastosowano metodę Powella - algorytm optymalizacji bezgradientowej wykorzystujący koncepcję kierunków sprzężonych. Metoda ta charakteryzuje się:
- Inicjalizacją z kierunkami bazowymi (wektory jednostkowe)
- Iteracyjną optymalizacją wzdłuż kolejnych kierunków
- Aktualizacją kierunków w celu uzyskania sprzężenia

**Parametry algorytmu:**
- Dokładność: $\varepsilon = 10^{-3}$
- Maksymalna liczba wywołań funkcji celu: $N_{max} = 10000$

### Metoda złotego podziału

Minimalizację jednowymiarową wzdłuż każdego kierunku przeprowadzono metodą złotego podziału, wykorzystującą własności złotej liczby $\phi = \frac{\sqrt{5} - 1}{2} \approx 0.618$.

### Zewnętrzna funkcja kary

Ograniczenia w problemie rzeczywistym uwzględniono stosując zewnętrzną funkcję kary z współczynnikiem kary $c = 10^6$:

$$\Phi(x) = f(x) + c \cdot \sum \max(0, g_j(x))^2$$

---

## Wyniki

### Walidacja implementacji

Przeprowadzono test walidacyjny dla belki o parametrach $l = 500$ mm, $d = 25$ mm:

| Parametr | Wartość obliczona | Wartość oczekiwana | Zgodność |
|----------|-------------------|-------------------|----------|
| Masa ($m$) | 2.189 kg | ~2.19 kg | ✓ |
| Ugięcie ($u$) | 36.217 mm | ~36.22 mm | ✓ |
| Naprężenie ($\sigma$) | 651.9 MPa | ~651.9 MPa | ✓ |

Uzyskane wartości potwierdzają poprawność implementacji równań fizycznych.

### Testowa funkcja celu

**Dla $w = 0$ (minimalizacja $f_2$):**
Wszystkie optymalizacje (niezależnie od wartości $a$) zbiegły do punktu $x^* \approx (-3, -3)$ z wartością $f_2^* \approx 0$, co odpowiada minimum globalnemu drugiego kryterium.

**Dla $w = 1$ (minimalizacja $f_1$):**
Wszystkie optymalizacje zbiegły do punktu $x^* \approx (3, 3)$ z wartością $f_1^* \approx 0$, co odpowiada minimum globalnemu pierwszego kryterium.

**Dla $w \in (0, 1)$:**
Rozwiązania tworzą front Pareto między dwoma ekstremalnymi punktami. Im większa wartość parametru $a$, tym silniejsza jest dominacja pierwszego kryterium, co wpływa na kształt frontu Pareto.

### Problem rzeczywisty - belka wspornikowa

**Kluczowe obserwacje:**

1. **Zbieżność do dolnego ograniczenia długości:** Większość rozwiązań zbiegła do $l^* = 200$ mm (minimalna dopuszczalna długość). Jest to zjawisko fizycznie poprawne, ponieważ:
   - Zarówno masa ($m \propto l$) jak i ugięcie ($u \propto l^3$) maleją wraz ze zmniejszaniem długości
   - Przy $l = 200$ mm i odpowiednim doborze średnicy wszystkie ograniczenia są spełnione
   - Minimalna długość reprezentuje rozwiązania Pareto-optymalne

2. **Wpływ współczynnika wagowego na średnicę:**
   - Dla $w = 0$ (minimalizacja ugięcia): $d^* \approx 48.7$ mm - duża średnica minimalizuje ugięcie
   - Dla $w = 1$ (minimalizacja masy): $d^* \approx 25.7$ mm - mniejsza średnica redukuje masę
   - Wartości pośrednie $w$ prowadzą do kompromisowych średnic

3. **Spełnienie ograniczeń:** Wszystkie znalezione rozwiązania spełniają ograniczenia:
   - Ugięcie: $u \leq 2.5$ mm ✓
   - Naprężenie: $\sigma \leq 300$ MPa ✓

---

## Dyskusja wyników

### Efektywność metody Powella

Metoda Powella okazała się bardzo efektywna w rozwiązywaniu problemów wielokryterialnych. Średnia liczba wywołań funkcji celu wynosiła:
- Problem testowy: ~200-300 wywołań na optymalizację
- Problem belki: ~400-450 wywołań na optymalizację

Metoda konsekwentnie znajdowała optima globalne lub punkty z frontu Pareto w zależności od wartości współczynnika wagowego.

### Metoda kryterium ważonego

Podejście z kryterium ważonym pozwoliło efektywnie wyznaczyć front Pareto poprzez 101 optymalizacji dla różnych wartości $w$. Uzyskane rozwiązania tworzą ciągły front reprezentujący wszystkie kompromisowe rozwiązania między dwoma kryteriami.

**Zalety metody:**
- Prosta implementacja - zamiana problemu wielokryterialnego na jednokryterialny
- Możliwość wykorzystania standardowych algorytmów optymalizacji
- Kontrola nad gęstością punktów na froncie Pareto (liczba wartości $w$)

**Ograniczenia:**
- Trudność w znalezieniu punktów na nieciągłych fragmentach frontu Pareto
- Konieczność wielokrotnego rozwiązywania problemu dla różnych $w$

### Wpływ parametru $a$ na front Pareto

Parametr $a$ znacząco wpływa na kształt frontu Pareto w problemie testowym:
- **$a = 1$:** Oba kryteria mają podobną wagę, front jest symetryczny
- **$a = 10$:** Pierwsze kryterium dominuje, front jest bardziej płaski
- **$a = 100$:** Silna dominacja pierwszego kryterium, front bardzo spłaszczony

### Zewnętrzna funkcja kary

Zastosowanie zewnętrznej funkcji kary z współczynnikiem $c = 10^6$ okazało się skuteczne w wymuszaniu spełnienia ograniczeń. Wszystkie znalezione rozwiązania dla problemu belki spełniały zarówno ograniczenia na ugięcie jak i naprężenie, jednocześnie optymalizując kryteria projektowe.

### Konsystencja punktów startowych

Zastosowanie identycznych punktów startowych dla różnych wartości parametru $a$ przy tym samym $w$ zapewniło uczciwą porównywalność wyników i pozwoliło na izolację wpływu parametru $a$ na kształt frontu Pareto.

---

## Wnioski

1. **Skuteczność metody Powella:** Metoda Powella w połączeniu z metodą złotego podziału okazała się bardzo skutecznym narzędziem do rozwiązywania problemów wielokryterialnych po zastosowaniu metody kryterium ważonego.

2. **Front Pareto:** Uzyskano kompletne fronty Pareto dla wszystkich rozpatrywanych przypadków, pokazujące zbiór optymalnych kompromisów między konkurującymi kryteriami.

3. **Problem belki:** Zbieżność rozwiązań do minimalnej dopuszczalnej długości jest fizycznie uzasadniona i reprezentuje optymalne rozwiązania przy danych ograniczeniach. Front Pareto pokazuje wymianę między masą a ugięciem poprzez zmianę średnicy belki.

4. **Zewnętrzna funkcja kary:** Metoda ta skutecznie wymusiła spełnienie wszystkich ograniczeń bez znaczącego zwiększenia złożoności obliczeniowej.

5. **Wpływ parametru $a$:** Parametr skali $a$ znacząco wpływa na kształt frontu Pareto, demonstrując jak różne wagi kryteriów mogą prowadzić do różnych zbiorów rozwiązań optymalnych.

6. **Praktyczne zastosowanie:** Wyniki pokazują, że optymalizacja wielokryterialna jest cennym narzędziem w problemach inżynierskich, gdzie konieczny jest kompromis między konkurującymi celami projektowymi.

---

## Załączniki

### Kod funkcji testowych

```cpp
// Pierwsze kryterium f1(x) = a*((x1-3)^2 + (x2-3)^2)
matrix ff5_f1(matrix x, matrix ud1, matrix ud2) {
    double a = m2d(ud1(0));
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    
    double result = a * (pow(x1 - 3.0, 2) + pow(x2 - 3.0, 2));
    return matrix(result);
}

// Drugie kryterium f2(x) = (1/a)*((x1+3)^2 + (x2+3)^2)
matrix ff5_f2(matrix x, matrix ud1, matrix ud2) {
    double a = m2d(ud1(0));
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    
    double result = (1.0 / a) * (pow(x1 + 3.0, 2) + pow(x2 + 3.0, 2));
    return matrix(result);
}

// Funkcja ważona: f(x) = w*f1(x) + (1-w)*f2(x)
matrix ff5T(matrix x, matrix ud1, matrix ud2) {
    double w = m2d(ud1(1));
    
    matrix f1 = ff5_f1(x, ud1, ud2);
    matrix f2 = ff5_f2(x, ud1, ud2);
    
    return w * f1 + (1.0 - w) * f2;
}
```

### Kod funkcji problemu rzeczywistego

```cpp
// Weighted objective for real problem: cantilever beam
// x(0) = diameter d [mm]
// x(1) = length l [mm]
// ud1(0) = weight w
// Returns weighted objective with penalties for constraint violations
matrix ff5R(matrix x, matrix ud1, matrix ud2) {
    double d_orig = m2d(x(0));  // diameter in mm
    double l_orig = m2d(x(1));  // length in mm
    double w = m2d(ud1(0)); // weight
    
    // Constants from K5.pdf
    const double P = 2000.0;        // force in N (2 kN)
    const double E = 120e9;         // Young's modulus in Pa (120 GPa)
    const double rho = 8920.0;      // density in kg/m^3
    const double u_max = 2.5;       // max deflection in mm
    const double sigma_max = 300e6; // max stress in Pa (300 MPa)
    
    // Bounds from K5.pdf
    const double d_min = 10.0;      // min diameter in mm
    const double d_max = 50.0;      // max diameter in mm
    const double l_min = 200.0;     // min length in mm
    const double l_max = 1000.0;    // max length in mm
    
    // External penalty coefficient
    const double penalty_coef = 1e6;
    
    // Penalty for bounds violations
    double f = 0.0;
    
    // Apply penalties for out-of-bounds, but clamp for calculation
    double d = d_orig;
    double l = l_orig;
    
    if (d_orig < d_min) {
        f += penalty_coef * pow(d_min - d_orig, 2);
        d = d_min;
    }
    if (d_orig > d_max) {
        f += penalty_coef * pow(d_orig - d_max, 2);
        d = d_max;
    }
    if (l_orig < l_min) {
        f += penalty_coef * pow(l_min - l_orig, 2);
        l = l_min;
    }
    if (l_orig > l_max) {
        f += penalty_coef * pow(l_orig - l_max, 2);
        l = l_max;
    }
    
    // Convert to SI units (meters)
    double d_m = d / 1000.0;
    double l_m = l / 1000.0;
    
    // Calculate mass (kg)
    double mass = rho * M_PI * pow(d_m / 2.0, 2) * l_m;
    
    // Calculate deflection (convert to mm)
    double u = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4));
    u = u * 1000.0;  // convert to mm
    
    // Calculate stress (Pa)
    double sigma = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
    
    // Criteria
    double f1 = mass;       // minimize mass
    double f2 = u;          // minimize deflection
    
    // Weighted objective
    f += w * f1 + (1.0 - w) * f2;
    
    // Add penalties for constraint violations (external penalty method)
    // Penalty for deflection constraint: u <= u_max
    if (u > u_max) {
        f += penalty_coef * pow(u - u_max, 2);
    }
    
    // Penalty for stress constraint: sigma <= sigma_max
    if (sigma > sigma_max) {
        f += penalty_coef * pow((sigma - sigma_max) / 1e6, 2);
    }
    
    return matrix(f);
}
```

### Funkcja lab5()

Główna funkcja implementująca pełne ćwiczenie lab5 znajduje się w pliku `main.cpp` (linie 778-955). Funkcja ta:
1. Przeprowadza 101 optymalizacji dla każdej wartości parametru $a$ w problemie testowym
2. Zapisuje wyniki do pliku `lab5_table1.csv` z zachowaniem identycznych punktów startowych dla tego samego $w$
3. Przeprowadza 101 optymalizacji dla problemu belki
4. Zapisuje wyniki do pliku `lab5_table2.csv`
5. Wyświetla postęp i podsumowanie obliczeń

**Wyniki optymalizacji** zostały przygotowane w formatach:
- `lab5_table1.csv` - Tabela 1: wyniki dla funkcji testowych (303 wiersze danych)
- `lab5_table2.csv` - Tabela 2: wyniki dla problemu belki (101 wierszy danych)

---

**Data wykonania:** [Data do uzupełnienia]

**Podpisy:**
- Jakub Sputo: ________________
- Jakub Stanula-Kaczka: ________________
- Kacper Strzesak: ________________
