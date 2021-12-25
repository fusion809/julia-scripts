#!/usr/bin/env julia
# For this recipe https://www.taste.com.au/recipes/vegetarian-shepherds-pie-cheesy-mash/tkk7u286?r=quickeasy&c=azmdgz4j/Easy%20vegetarian%20dinners&h=Quick%20%20easy
cats = ["Energy kJ";
	"Protein (g)";
	"Fat, total (g)";
	"Fat, saturated (g)";
	"Carbohydrate (g)";
	"Sugars (g)";
	"Dietary fibre (g)";
	"Sodium (mg)"]; # Per 100 g
onionPer100g = [127;
		1.7;
		0.1;
		0;
		4.6;
		4.6;
		2.1;
		11] # Per Woolies website
onionPer100gNUTTAB = [139;
		      1.3;
		      0;
		      0;
		      5.8;
		      5.8;
		      2.7;
		      8];
onionMass = 1.42 # Per NUTTAB
celeryPer100g = [62;
		 0.6;
		 0.1;
		 0;
		 1.2;
		 1.2;
		 1.4;
		 98];
celeryMass = 2/7 * 3; # Estimated as the 300g celery sticks punnet image looks like 7 sticks at 300 g, and we need 2
carrotPer100g = [132;
		 0.8;
		 0.1;
		 0;
		 5;
		 5;
		 3.9;
		 40]
carrotMass = 1.29;
zucchiniMass = 2* 1.95;
zucchiniPer100g = [61;
		   0.8;
		   0.3;
		   0;
		   1.6;
		   1.6;
		   1.2;
		   1];
kentPumpkinMass = 5;
kentPumpkinPer100g = [186;
		      1.4;
		      0.2;
		      0;
		      8.1;
		      6.1;
		      2.6;
		      1];
finelyChoppedTomatoesMass = 2*4;
finelyChoppedTomatoesPer100g = [110;
				1.2;
				1;
				0.5; # The fat and sat fat content are guesses, as <1g is given for both
				3.9;
				2.8;
				1.2; # Guess also, based on the 1.2 g per 100 g for raw tomatoes
				120];
veggieStockPer100g = [12.4;
	       0.15;
	       0;
	       0;
	       0.53;
	       0.36;
	       0.3;
	       329];
veggieStockMass = 4.375;
redLentilsMass = 1.25;
redLentilsPer100g = [1390;
		     25.4;
		     2.4;
		     1; # Estimate as it says <1 g
		     43.5;
		     1.2;
		     15.9;
		     5] # Estimate as it says <5 mg
soySaucePer100g = [113;
		   1; # Estimate as its says <1g
		   0;
		   0;
		   7;
		   0;
		   0; # Estimate as it's omitted
		   7300]; # For https://www.woolworths.com.au/shop/productdetails/841770/silver-swan-soy-sauce
soySauceMass = 0.18;
brownLentilPer100g = [373;
		      6.8;
		      0.3;
		      0.1;
		      12.6;
		      0.2;
		      6.4;
		      274];
brownLentilMass = 2 * 1.2;
whitePotatoMass = 13;
whitePotatoPer100g = [281;
		      2.3;
		      0.1;
		      0; # Estimated as said <0.1
		      12.9;
		      0.6;
		      1.8;
		      6];
appleCiderVinegarPer100g = [65;
			 0.1;
			 0; # All these zero fields are omitted
			 0;
			 0;
			 0;
			 0;
			 5];
appleCiderVinegarMass = 0.18;
garlicCloveMass = 0.08; # DuckDuckGo gave 4g as the mass of a clove and we're using 2 of them
garlicClovePer100g = [522;
		      6.1;
		      2.8;
		      0.35*2.8;
		      10.2;
		      1.5;
		      1.6;
		      8];
milkPer100g = [292;
	       3.5;
	       3.4;
	       2.3;
	       6.3;
	       6.3;
	       0;
	       37];
milkMass = 1.25;
butterMass = 0.8;
butterPer100g = [3050;
		 0.5; # Estimate as it says <1g;
		 82.1;
		 53.1;
		 0.5; # Also an estimate as says <1g
		 0.5; # Also an estimate as says <1g
		 0;
		 550];
cheeseMass = 1.25;
# Going to use Mozzarella cheese, despite the recipe calling for tasty, as it's healthier
cheesePer100g = [1170;
		25.9;
		19.2;
		12.7;
		0.5; # <1g is listed
		0.5; # <1g is listed
		0;
		483];
totalNutr = cheesePer100g * cheeseMass + butterPer100g * butterMass + milkMass * milkPer100g + garlicClovePer100g * garlicCloveMass + appleCiderVinegarMass * appleCiderVinegarPer100g + whitePotatoPer100g * whitePotatoMass + brownLentilMass * brownLentilPer100g + soySauceMass * soySaucePer100g + redLentilsPer100g * redLentilsMass + veggieStockMass * veggieStockPer100g + finelyChoppedTomatoesMass * finelyChoppedTomatoesPer100g + kentPumpkinMass * kentPumpkinPer100g + zucchiniMass * zucchiniPer100g + carrotMass * carrotPer100g + celeryMass * celeryPer100g + onionPer100g * onionMass
perServeNutr = totalNutr / 6
println("Per serve the vegetarian shepherd's pie has: ")
println("Energy (kJ): ", perServeNutr[1])
println("Protein (g): ", perServeNutr[2])
println("Fat (g): ", perServeNutr[3])
println("Saturated fat (g): ", perServeNutr[4])
println("Carbohydrate (g): ", perServeNutr[5])
println("Sugars (g): ", perServeNutr[6])
println("Dietary fibre (g): ", perServeNutr[7])
println("Sodium (mg): ", perServeNutr[8])