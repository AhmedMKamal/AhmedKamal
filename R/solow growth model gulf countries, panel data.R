#first install packages

install.packages('WDI'); install.packages('plm')
install.packages("stargazer");install.packages('tidyverse');install.packages('psych')


# Set the indicator codes for GDP per capita, population growth, and gross capital formation

library(WDI);library(stargazer)
library(tidyverse);library(psych)


gdp_capita <- "NY.GDP.PCAP.KD"
labor_force_code <- "SL.TLF.TOTL.IN"
gc_form_gdp_code <- "NE.GDI.TOTL.ZS"
gdp_code <- "NY.GDP.MKTP.KD"

# Set the country codes for Arab Gulf States
gulf_states <- c("BHR", "KWT", "OMN", "QAT", "SAU", "ARE", "YEM")

# Set the start and end years for the data
start_year <- 2000
end_year <- 2022

# Download the data using the WDI package
gdp_capita <- WDI(country = gulf_states, indicator = gdp_capita, start = start_year, end = end_year)
labor_force <- WDI(country = gulf_states, indicator = labor_force_code, start = start_year, end = end_year)
capital_of_gdp <- WDI(country = gulf_states, indicator = gc_form_gdp_code, start = start_year, end = end_year)
gdp_data <- WDI(country = gulf_states, indicator = gdp_code, start = start_year, end = end_year)

# Clean the data using dplyr
gdp_capita <- gdp_capita %>% 
  select(country, year, NY.GDP.PCAP.KD) %>% 
  rename(gdp_per_capita = NY.GDP.PCAP.KD)

labor_force <- labor_force %>% 
  select(country, year, SL.TLF.TOTL.IN) %>% 
  rename(labor_force = SL.TLF.TOTL.IN)

capital_of_gdp <- capital_of_gdp %>% 
  select(country, year, NE.GDI.TOTL.ZS) %>% 
  rename( capital_of_gdp = NE.GDI.TOTL.ZS)

gdp_data <- gdp_data %>% 
  select(country, year, NY.GDP.MKTP.KD) %>% 
  rename(gdp_constant_usd = NY.GDP.MKTP.KD)

# Merge the data into one data frame
gulf_data <- merge(gdp_data, labor_force, by = c("country", "year"))
gulf_data <- merge(gulf_data, capital_of_gdp, by = c("country", "year"))
gulf_data<- merge(gulf_data, gdp_capita, by = c("country", "year"))

#visual representation

gulf_data <- gulf_data%>%
  group_by( country )%>%
  mutate(LFGR = 100* (labor_force - lag(labor_force ) )/
           lag(labor_force ) )%>%
  view()


# Create plots and charts
ggplot(gulf_data, aes(x = log(capital_of_gdp), y = log(gdp_per_capita), color=country )  ) +
  geom_point() +
  labs(x = "log Capital Formation", y = "log GDP Per Capita") +
  theme_classic()+ geom_smooth(method = "lm", se = FALSE)

ggplot(gulf_data, aes(x = log(LFGR), y = log(gdp_per_capita), color = country)) +
  geom_point() +
  labs(x = "Labor Force Growth Rate ", y = "llog GDP Per Capita", color = "Country") +
  theme_classic()+ geom_smooth(method = "lm", se = FALSE)

ggplot(gulf_data, aes(x = year, y = gdp_per_capita, color = country)) +
  geom_point() +
  labs(x = "Year", y = "GDP Per capita", color = "Country") +
  geom_line()+
  theme_classic()


#prepare the data for regression analysis

#remove Yemen, Rep from our data
library(plm)
gulf_data_subset <- subset(gulf_data, country != "Yemen, Rep.")
gulf_data_subset<- pdata.frame(gulf_data_subset, index =c("country", "year") )
head(gulf_data_subset)
min(gulf_data_subset$LFGR, na.rm = T)

#estimate models
Pooled<- plm( log(gdp_per_capita) ~ log(capital_of_gdp) + log(LFGR+10), data=gulf_data_subset,
              model = 'pooling' ) 

Fixed<- plm( log(gdp_per_capita) ~ log(capital_of_gdp) + log(LFGR+10), data=gulf_data_subset,
             model = 'within' ) 

Random<- plm( log(gdp_per_capita) ~ log(capital_of_gdp) + log(LFGR+10), data=gulf_data_subset,
              model = 'random' ) 


GMM<-  pgmm(log(gdp_per_capita) ~ lag( log(gdp_per_capita) ) + log(capital_of_gdp) + log(LFGR+10)
            | lag( log(capital_of_gdp),2:6 ) + lag( log(LFGR+10),2:6 ),
            effect = "twoways",
            data = gulf_data_subset, 
            mdeaths="twosteps", transformation = "ld" , collapse = F)

summary(Pooled)
summary(Fixed)
summary(Random)
summary(GMM)
#hausman test 
phtest(Fixed, Random)


#Organize the outputs into a table
library(stargazer)
stargazer(Pooled, Fixed, Random, GMM, type = "text",
          column.labels = c("Pooled", "Fixed","Random","GMM") )


#descriptive statistics
describe(cbind(Lnn=log(gulf_data_subset$LFGR +10 ),
               LnGDP=log(gulf_data_subset$gdp_constant_usd),
               Lns=log(gulf_data_subset$capital_of_gdp)    ) )







