
# include("/Users/tarikz/Documents/Code/01-julia/05-PlotData.jl")


using CSV

#data = CSV.read("/Users/tarikz/Documents/Code/01-julia/05-NAEXKP01TRQ657S.csv")
data_lux = CSV.read("/Users/tarikz/Documents/Code/01-julia/05-b3019.csv", datarow = 2, delim = ',')
# Use 'showcols(data_lux)' to see the variable names


### Method 1: Using Plots
#using Plots
#plotly()
#plot( data.DATE, data.NAEXKP01TRQ657S, linewidth=2, title="GDP Turkey"  )

### Method 2: Using Gadfly (visit: http://gadflyjl.org)
using Gadfly
#draw(SVG("/Users/tarikz/Documents/Code/01-julia/output.svg", 6inch, 3inch), plot(data.DATE, data.NAEXKP01TRQ657S))
#plot(data.DATE, data.NAEXKP01TRQ657S)

### Reorganise the data (STATEC organises time series from latest to oldest)

#│ Row │ variable                                              │ eltype   │ nmissing │ first          │ last         │
#│     │ Symbol                                                │ DataType │ Int64    │ Any            │ Any          │
#├─────┼───────────────────────────────────────────────────────┼──────────┼──────────┼────────────────┼──────────────┤
#│ 1   │ Spécification                                         │ String   │ 0        │ septembre 2018 │ janvier 1995 │
#│ 2   │ 1. Salariés résidents sortants - données cvs          │ Float64  │ 0        │ 12551.5        │ 8628.11      │
#│ 3   │ 2. Salariés frontaliers entrants - données cvs        │ Float64  │ 0        │ 1.93647e5      │ 53773.7      │
#│ 4   │ 3. Emploi salarié national - données cvs              │ Float64  │ 0        │ 2.43207e5      │ 1.49295e5    │
#│ 5   │ 4. Emploi salarié intérieur - données cvs (3 + 2 - 1) │ Float64  │ 0        │ 4.24302e5      │ 1.9444e5     │
#│ 6   │ 5. Emploi non-salarié national - données cvs          │ Float64  │ 0        │ 22144.2        │ 16912.3      │
#│ 7   │ 6. Emploi non-salarié intérieur - données cvs         │ Float64  │ 0        │ 27091.1        │ 17462.6      │
#│ 8   │ 7. Emploi total national - données cvs (3 + 5)        │ Float64  │ 0        │ 2.65351e5      │ 1.66207e5    │
#│ 9   │ 8. Emploi total intérieur - données cvs (4 + 6)       │ Float64  │ 0        │ 4.51394e5      │ 2.11903e5    │
#│ 10  │ 9. Nombre de chômeurs - données cvs                   │ Float64  │ 0        │ 15027.0        │ 4331.61      │
#│ 11  │ 10. Population active - données cvs (7 + 9)           │ Float64  │ 0        │ 2.80378e5      │ 1.70538e5    │
#│ 12  │ 11. Taux de chômage (en %) - données cvs (9 / 10)     │ Float64  │ 0        │ 5.35954        │ 2.53996      │



nrows, ncols = size(data_lux)
df = zeros(nrows, 1)
for r in 0:nrows-1
     df[r+1, 1] = data_lux[nrows-r, 12]   
end
# Working with dates: https://en.wikibooks.org/wiki/Introducing_Julia/Working_with_dates_and_times
dr = Dates.Date(1995,1):Dates.Month(1):Dates.Date(2018,9)

using Gadfly, Compose
# Visit: http://gadflyjl.org
plot(dr, df, label="empl")


