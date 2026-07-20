#!/usr/bin/env python3
# Viraj Arora

import argparse
import math
import os
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore
import scipy as sp # type: ignore

class findPlasmidLossRate:


    def __init__(self, od_pcb, viability, prop_dark, treatment='Carbenicillin'):

        self.od_pcb = od_pcb
        self.viability = viability
        self.prop_dark = prop_dark
        self.treatment = treatment

        #loads dataframes from the input files
        self.odpcb_df   = pd.read_csv(self.od_pcb, sep='\t')
        self.death_df   = pd.read_csv(self.viability, sep='\t')

        # gets the average of the replicates for each treatment in each dataframe and stores it as a new dataframe
        self.avg_odpcb_df = self.averageReplicates(self.odpcb_df)
        self.avg_via_df = self.averageReplicates(self.death_df)

        # prop_dark can be a filepath to a csv, a numeric proportion, or a list of (prop, treatment) tuples
        if isinstance(self.prop_dark, list):
            carb = [p for p, t in self.prop_dark if t == 'Carbenicillin']
            kan  = [p for p, t in self.prop_dark if t == 'Kanamycin']
            n = max(len(carb), len(kan))
            carb = carb + [np.nan] * (n - len(carb))
            kan  = kan  + [np.nan] * (n - len(kan))
            self.proportions_df = pd.DataFrame(
                {'Carbenicillin': carb, 'Kanamycin': kan},
                index=[f'delta{i+1}(avg)' for i in range(n)]
            )
        elif isinstance(self.prop_dark, (int, float)):
            self.proportions_df = pd.DataFrame(
                {self.treatment: [float(self.prop_dark)]},
                index=['delta1(avg)']
            )
        else:
            self.propdark_df = pd.read_csv(self.prop_dark, index_col=0) # prop_dark file is a csv where all relevant data is in rows labelled "delta 6 #A$"
            self.proportions_df = self.cleanPropDark(self.propdark_df) # uses helper method to pull necessary data from prop_dark file, and averages over replicates

        #calculates the doubling time to be used for all calculations of plasmid loss rate 
        self.d_time = self.getDoublingTime()
    
    def cleanPropDark(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Takes the experimental data of rows labeled with 'delta #X' and returns a dataframe with final plasmid loss proportions for each antibiotic treatment. Colony counts are averaged over each replicate, denoted by shared number. 
        '''

        df_new = df[df.index.str.contains(r'delta\s\d', regex=True, na=False)]

        count = sum(df_new.index.str.contains(r'delta\s\d', regex=True))

        length_new = int(count / 5)

        cleaned_data = pd.DataFrame({'Carbenicillin': [0] * length_new, 'Kanamycin' : [0] * length_new}, 
                                    index=[f'row{i}' for i in range(length_new)]
                        )

        # converts to numeric values(1) -> sum plating counts over each replicate and calculate overall proportion (2) ->  get average proportions for replicates of each treatment(3)

        for col in df_new.columns: #(1)
            df_new[col] = pd.to_numeric(df_new[col], errors='coerce') 

        # this step is done with respect to the exact dimensions of the data given, so needs to be modified if data format changes
        for ind in range(len(df_new)): #(3)
            if ind % 5 == 4: # 5 works due to it being th enumber of trials per replicate
                count = int(ind / 5)
                temp_row = df_new.iloc[ind - 4:ind + 1].sum() #(2)

                prev_val = None

                # implementation specific to format given, just gets dark colonies and total colonies for each treatment
                for ind, val in enumerate(temp_row.values):
                    if ind % 2 == 0:
                        prev_val = val
                    elif ind % 2 == 1:
                        cleaned_data.iloc[int(count), (ind-1) // 2] = prev_val / val

                cleaned_data = cleaned_data.rename(index={cleaned_data.index[count]: 'delta' + str(count + 1) + '(avg)'})
                
        return cleaned_data
        
        

    def averageReplicates(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Averages replicate columns in a dataframe by shared treatment tag.
        '''
        # turns empty cells become NaN (avoids errors in averaging)
        df = df.replace(r'^\s*$', np.nan, regex=True)
        df.columns = df.columns.str.replace(r'\.\d+$', '', regex=True)

        # Group and average replicate columns and average them
        grouped = df.groupby(df.columns, axis=1).mean()
        return grouped

    def getDoublingTime(self):
        '''Gets the average doubling time of the population from OD data - uses the "Top10 / pGFPck Δ6" strain grown under carbenicillin as the reference culture'''
        df = self.avg_odpcb_df.copy()
        col = "Top10 / pGFPck Δ6"

        time_col = df.columns[0]

        # Calculates doubling time from the maximum and minimum OD values from the growth data and the time between them, as growth is monotonic. 
        # note : this can be changed at a later point in prefereence of a more logical interval selection rule
        max_time = df[time_col].max()
        od_fchange = df[col].max() / df[col].min()
        dtime = max_time / (math.log(od_fchange) / math.log(2))

        return dtime
    
    def getLosses(self):
        '''Get plasmid loss rates for each replicate and treatment'''
        loss_rates = {'Kanamycin': [], 'Carbenicillin': []} # empty dictionary to store loss rates

        for ind, row in self.proportions_df.iterrows():
            if ind != 'total_avg':
                rates_curr = self.getLoss(row)
                for key in rates_curr:
                    if key.split('.')[0] == 'Carbenicillin':
                        loss_rates['Carbenicillin'].append(rates_curr[key])
                    elif key.split('.')[0] == 'Kanamycin':
                        loss_rates['Kanamycin'].append(rates_curr[key])

        
        return loss_rates

    def getLoss(self, row: pd.Series) -> dict: 
        'compute the loss rate for one replicate of each treatment'

        beta = math.log(2) / self.d_time
        rates = {}
        for treatment in row.index:
            #strips .x from treatment labels in propdark_df to match them with the labels in the viability dataframe
            treat = treatment.split('.')[0]
            prop_dark = row[treatment]
            if pd.isna(prop_dark):
                continue
            # calculates C / F, assumes data is in form of percentage of dark cells
            prop = float(prop_dark)
            R = (1 - prop) / prop

            # gets a linear interpolation of the survival curve for current treatment label 
            time_col = self.avg_via_df['Hours']
            via_prop = self.avg_via_df[treat]
            valid_prop = via_prop.notna()
            valid_time = time_col[valid_prop]
            valid_prop = via_prop[valid_prop]

            T = valid_time.max() # gets the length of the survival curve
            p = lambda t: np.interp(t, valid_time, valid_prop) # defines the survival function as a linear interpolation of the viability data for the current treatment

            # calculates the loss rate using the formula derived in the notebook, wiht scipy's quad and brentq (root finder) functions
            # note: the step subtracts by the final plasmid loss ratio observed at plating to turn the numerical method into a root finding problem (i.e. solve for L such that the equation equals 0)
            rates[treatment] = sp.optimize.brentq(lambda L: 
                np.exp((beta-L)*T) / (L * sp.integrate.quad(lambda t: np.exp((beta-L)*t) * p(T-t), 0, T)[0]) - R, 
                1e-10, 10)
            rates[treatment] = rates[treatment] * self.d_time # converts loss rate to per division by multiplying by doubling time

        return rates

    def plotLosses(self, rates):
        '''Bar plot of loss rates by treatment. Call only when a graph is requested.'''
        plt.bar(rates.keys(), [v[0] for v in rates.values()], yerr=[v[1] for v in rates.values()], capsize=5, color='gray', edgecolor='black')
        plt.ylabel("Loss Rate (per division)")
        plt.xlabel("Treatment")
        plt.title("Plasmid Loss Rates by Antibiotic")
        plt.show()





def parseProportion(value: str):
    '''Parse the prop_dark CLI argument. Returns (prop, treatment) where prop is a float or a filename string.
    A trailing "k" (case-insensitive) marks the value as Kanamycin data; otherwise Carbenicillin is assumed.
    Whitespace is ignored throughout.'''
    stripped = ''.join(value.split())
    treatment = 'Carbenicillin'
    if stripped and stripped[-1].lower() == 'k':
        treatment = 'Kanamycin'
        stripped = stripped[:-1]
    # ratio form "x/y"
    if '/' in stripped:
        num, den = stripped.split('/', 1)
        return float(num) / float(den), treatment
    # decimal form
    try:
        return float(stripped), treatment
    except ValueError:
        # otherwise treat as a filename (use the original value so filenames with spaces/dots are preserved)
        return value, treatment


def main():
    parser = argparse.ArgumentParser(description="Compute plasmid loss rates from one or more proportions, or from a colony count file.")
    parser.add_argument("prop_dark", nargs='+', help="One or more proportions of dark colonies (decimal like 0.12, or ratio like 3/25; append 'k' for Kanamycin, default is Carbenicillin), or a single csv filename in this folder.")
    parser.add_argument("-g", "--graph", action="store_true", help="Display a bar plot of the loss rates.")
    parser.add_argument("-o", "--output", nargs='?', const='results', default=None, help="For proportion inputs: write results to a .txt file in the script folder. Optional name follows -o (default 'results'); a trailing .txt is ignored. Ignored in csv file mode, which always writes <filename>_results.txt.")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    odpcb_file = os.path.join(script_dir, "OD_values_growth_curve_uv_vs_ck_070126.txt") # fixed file to load growth rate data for plasmid-containing bacteria
    viability_file = os.path.join(script_dir, "NORMALIZED_CFU_vs_time.txt") # fixed file to load viability data for plasmid-free bacteria in order to estimate survival rates of bacteria who lose their plasmid

    # If a single arg is passed and it parses as a filename, use file mode. Otherwise treat all args as proportions.
    parsed = [parseProportion(v) for v in args.prop_dark]
    proportion_mode = not (len(parsed) == 1 and isinstance(parsed[0][0], str))
    if not proportion_mode:
        prop_dark = os.path.join(script_dir, parsed[0][0])
        treatment = parsed[0][1]
    else:
        prop_dark = [(float(p), t) for p, t in parsed]
        treatment = 'Carbenicillin' # unused when prop_dark is a list

    plasmid_loss_obj = findPlasmidLossRate(odpcb_file, viability_file, prop_dark, treatment=treatment)

    rates = plasmid_loss_obj.getLosses()
    # handles data for each treatment to get the averages and standard deviation accordingly
    for head in list(rates.keys()):
        if len(rates[head]) == 0:
            del rates[head]
        else:
            rates[head] = [np.mean(rates[head]), np.std(rates[head], ddof = 0)]

    lines = []
    if proportion_mode:
        lines.append("Inputs:")
        by_treatment = {}
        for raw, (_, t) in zip(args.prop_dark, parsed):
            by_treatment.setdefault(t, []).append(raw)
        for t, raws in by_treatment.items():
            lines.append(f"  {t}: {', '.join(raws)}")

    lines.append("Plasmid Loss Rates (per division):")
    for head in rates:
        if proportion_mode:
            lines.append(f"{head}: {rates[head][0]:.5f}")
        else:
            lines.append(f"{head}: {rates[head][0]:.5f} ± {rates[head][1]:.5f}")

    output_text = "\n".join(lines)
    print(output_text)

    # File output: csv mode always writes; proportion mode writes only when -o is given.
    out_name = None
    if not proportion_mode:
        stem = os.path.splitext(os.path.basename(parsed[0][0]))[0]
        out_name = f"{stem}_results"
    elif args.output is not None:
        out_name = args.output
        if out_name.lower().endswith('.txt'):
            out_name = out_name[:-4]

    if out_name is not None:
        out_path = os.path.join(script_dir, out_name + '.txt')
        with open(out_path, 'w') as f:
            f.write(output_text + "\n")

    if args.graph:
        plasmid_loss_obj.plotLosses(rates)


if __name__ == "__main__":
    main()