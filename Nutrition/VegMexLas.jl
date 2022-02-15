# This is for the recipe https://www.taste.com.au/recipes/vegetarian-mexican-lasagne/569f18fb-3163-4b92-9ff5-3f824e7ad3a5
# Except we've substituted Mozzarella cheese for tasty cheese as it's healthier
# Nutrition info order:
# - Energy (kJ)
# - Protein (g)
# - Fat (g)
# - Saturated fat (g)
# - Carbohydrates (g)
# - Sugars (g)
# - Dietary fibre (g)
# - Sodium (mg)
# URL: https://www.woolworths.com.au/shop/productdetails/780390/perfect-italiano-grated-mozzarella-cheese
mozCheesePer100g = [1250,
26.1,
22.0,
13.2,
0.1,
0.1,
0,
527]
mozCheeseServ = 4/4 * 2.5/4 # Per serving size in units of 100g

# Multigrain tortilla
# URL: https://www.woolworths.com.au/shop/productdetails/333917/mission-multigrain-tortillas
multGrainTortPer100g = [1280,
8.6,
9.1,
4.3,
44.5,
5.2,
4.4,
1070]
multGrainServ = 5*0.48/4 # Tortillas are 48 g each and we're using 5 giving us 4 serves

# Edgell Corn Kernels 300 g (closest available at time of writing to 310g)
# URL: https://www.woolworths.com.au/shop/productdetails/966533/edgell-corn-kernels
cornPer100g = [261,
2.3,
1.6,
0.3,
8.0,
4.0,
3.3,
251]
cornServ = 0.95*2.7/4 # Some of the can is just liquid, not the actual corn.

# Woolworths Red Kidney Beans No Added Salt 420G
# URL: https://www.woolworths.com.au/shop/productdetails/259401/woolworths-red-kidney-beans-no-added-salt
kidneyBeanPer100g = [423,
7.8,
0.5,
0.1,
12.9,
0.6,
6.5,
5]
kidneyBeanServ = 0.75*3/4

# Old El Paso Mexican Taco Spice Mix 30G
# URL: https://www.woolworths.com.au/shop/productdetails/93444/old-el-paso-mexican-taco-spice-mix
tacoSeasonPer100g = [1160,
7.5,
3.3,
0.4,
50.3,
33.2,
6.7,
6620]
tacoSeasonServ = 0.3/4

# Old El Paso Mexican Reduced Salt Taco Spice Mix 30G
# URL: https://www.woolworths.com.au/shop/productdetails/65083/old-el-paso-mexican-reduced-salt-taco-spice-mix
tacoSeasonRdSaltPer100g = [1250,
7.5,
3.2,
0.4,
55.7,
33.3,
6.8,
3990]
tacoSeasonRdSaltServ = 0.3/4

# Dolmio Extra Tomato, Onion & Garlic Salt Reduced Pasta Sauce 500g
# URL: https://www.woolworths.com.au/shop/productdetails/758316/dolmio-extra-tomato-onion-garlic-salt-reduced-pasta-sauce
pastaSaucePer100g = [190,
1.8,
0.1,
0.1,
9,
4.9,
2.3,
159]
pastaSauceServ = 5/4

# Onion Brown Each
# URL: https://www.woolworths.com.au/shop/productdetails/144329/onion-brown
onionPer100g = [127,
1.7,
0.1,
0,
4.6,
4.6,
2.1,
11]
# Serving size approximated from 1 brown onion called for in recipe and NUTTAB indicated the size is about 142g
onionServ = 1.42/4

garlicCloveServ = 0.08/4; # DuckDuckGo gave 4g as the mass of a clove and we're using 2 of them
# Nutritional info from NUTTAB https://www.foodstandards.gov.au/science/monitoringnutrients/afcd/Pages/fooddetails.aspx?PFKID=F004193
garlicClovePer100g = [522,
6.1,
2.8,
0.35*2.8,
10.2,
1.5,
16.9,
8];

# Green Zucchini Each
# URL: https://www.woolworths.com.au/shop/productdetails/170225/green-zucchini
zucchiniServ = 1.95/4; # NUTTAB estimated size as 195 g for green zucchini
zucchiniPer100g = [61,
0.8,
0.3,
0,
1.6,
1.6,
1.2,
1];

# Red Capsicum https://www.woolworths.com.au/shop/productdetails/135306/red-capsicum
redCapsicumPer100g = [106,
1.5,
0.2,
0,
3.5,
3.5,
1.8,
2];
redCapsicumServ = 2.44/4; # NUTTAB estimates red capsicum size is 244 g

nutrServ = mozCheesePer100g * mozCheeseServ + multGrainTortPer100g * multGrainServ + cornPer100g * cornServ + kidneyBeanPer100g * kidneyBeanServ + tacoSeasonRdSaltPer100g * tacoSeasonRdSaltServ + pastaSaucePer100g * pastaSauceServ + onionPer100g * onionServ + garlicClovePer100g * garlicCloveServ + zucchiniPer100g * zucchiniServ + redCapsicumPer100g * redCapsicumServ

# From https://stackoverflow.com/a/46292069/1876983
# With some modification
function signifChop(num, digits)
    if num == 0.0 then
        return num
    else
        e = ceil(log10(abs(num)))
        scale = 10^(digits - e)
        res = trunc(num * scale) / scale
        if res == round(res)
            Int64(res) # Show as integer without .0, if integer
        else
            return res
        end
    end
end

# Print food-specific info for specific nutrient
function printNutr(nutr, no, unit)
    str = nutr * " per serve " * (" "^(13-length(nutr))) * "= "
    unit = (" "^ (2 - length(unit))) * unit
    sigfigs = 4
    val = string(signifChop(nutrServ[no], sigfigs))
    lenVal = length(val)
    val = val * (" "^(sigfigs+1-lenVal))

    println(str, val, " " * unit)
end

printNutr("Energy", 1, "kJ")
printNutr("Protein", 2, "g")
printNutr("Fat", 3, "g")
printNutr("Saturated fat", 4, "g")
printNutr("Carbohydrate", 5, "g")
printNutr("Sugars", 6, "g")
printNutr("Dietary fibre", 7, "g")
printNutr("Sodium", 8, "mg")