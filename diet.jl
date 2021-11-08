# Weet-bix
weetbixEnergy = 492 * 2; # Four bricks
# https://www.woolworths.com.au/shop/productdetails/200695/sanitarium-weet-bix-breakfast-cereal⎈

soyMilkEnergy = 618 * 4; # 1 litre of vanilla soymilk
# https://www.woolworths.com.au/shop/productdetails/5535/sanitarium-so-good-long-life-vanilla-bliss-soy-milk

bananaEnergy = 386 * 2; # 2 bananas
# https://www.foodstandards.gov.au/science/monitoringnutrients/afcd/Pages/fooddetails.aspx?PFKID=F000262

miloEnergy = 310 * 2; # At most I might use two servings
# plant-based milo https://www.woolworths.com.au/shop/productdetails/94324/nestle-milo-plant-based

carrotEnergy = 190 * 2; # Two carrots
# 129 g carrot https://www.foodstandards.gov.au/science/monitoringnutrients/afcd/Pages/fooddetails.aspx?PFKID=F002276

grapeEnergy = 1610; # 500g of grapes
# https://www.foodstandards.gov.au/science/monitoringnutrients/afcd/Pages/fooddetails.aspx?PFKID=F004259

# soy and linseed bread
soyLinseedEnergy = 1070*2 ; # Serving 

# Baked beans
beansEnergy = 444*3; # Can of baked beans
# https://www.woolworths.com.au/shop/productdetails/47682/woolworths-baked-beans-in-tomato-sauce

totalEnergy = weetbixEnergy + soyMilkEnergy + bananaEnergy + miloEnergy + carrotEnergy + grapeEnergy + soyLinseedEnergy + beansEnergy;

# Energy requirements calculator estimates this energy requirement
# https://www.eatforhealth.gov.au/webform/daily-energy-requirements-calculator⎈
reqEnergyLower = 9787;
reqEnergyHigher = 11186;

diffLower = reqEnergyLower - totalEnergy;
diffHigher = reqEnergyHigher - totalEnergy;
println("You have consumed at least ", diffLower, "kJ less energy than you need");
println("You have consumed at most ", diffHigher, "kJ less energy than you need");

