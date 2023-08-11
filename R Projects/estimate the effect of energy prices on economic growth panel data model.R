#install packages

install.packages('psych')
install.packages('plm')
install.packages('readxl')
install.packages('tidyverse')
install.packages('stargazer')

library(readxl);library(plm);library(psych)
library(stargazer); library(tidyverse)

#Import Data

data <- read_excel("D:/Economitrics/Eurropean data.xlsx", 
                   sheet = "from2007")

#transform the data to panal data frame

pdata<- pdata.frame(data, index = c("country","years"))

#descriptive statistics

describe( log( cbind(select(pdata,GDP:TR,S ) ,pdata$LF+0.1+0.05) ) )

#visual representation
ggplot(data, mapping = aes(x=log(COP), y=log(GDP) ) )+
  geom_point(aes(color=country) )+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

ggplot(data, mapping = aes(x=log(NGP) , y=log(GDP) ) )+
  geom_point(aes(color=country) )+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

ggplot(data, mapping = aes(x=log(FEC) , y=log(GDP) ) )+
  geom_point(aes(color=country) )+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

ggplot(data, mapping = aes(x=log(TR) , y=log(GDP) ) )+
  geom_point(aes(color=country) )+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

ggplot(data, mapping = aes(x=log(S) , y=log(GDP) ) )+
  geom_point(aes(color=country) )+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

ggplot(data, mapping = aes(x=log(LF+0.05+0.1) , y=log(GDP) ) )+
  geom_point(aes(color=country) )+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

ggplot(data, mapping = aes(x=years, y=GDP)) + 
  geom_line() + facet_wrap(country~.)+ 
  theme_minimal()+ scale_x_continuous(labels = NULL)


#estimate the panel data models

Pooled<- plm(log(GDP)~log(COP)+log(NGP)+log(FEC)+log(S)+log(TR)+log(LF+0.1+0.05),
             data = pdata ,model = "pooling")

Fixed<- plm(log(GDP)~log(COP)+log(NGP)+log(FEC)+log(S)+log(TR)+log(LF+0.1+0.05),
            data = pdata ,model = "within")

Random<- plm(log(GDP)~log(COP)+log(NGP)+log(FEC)+log(S)+log(TR)+log(LF+0.1+0.05),
             data = pdata ,model = "random")

summary(plm(log(GDP)~log(COP)+log(NGP)+log(FEC)+log(S)+log(TR)+log(LF+0.1+0.05),
            data = pdata ,model = "within") )

GMM<- pgmm(log(GDP) ~ lag(log(GDP),2) + log(COP) + log(NGP) + log(FEC) +
             log(S) + log(TR) + log(LF+0.1+0.05) | lag(log(GDP), 3:99) +
             lag(log(COP), 2:99) + lag(log(NGP), 2:99) + lag(log(FEC), 2:99) +
             lag(log(S), 2:99) + lag(log(TR), 2:99) + lag(log(LF+0.1+0.05), 2:99),
           data = pdata, effect = "individual", model = "onestep", transformation = "d", collapse = T)

#Hausman Test
phtest(Fixed,Random)

#Organize the outputs into a table
stargazer(list(Pooled,Fixed,Random,GMM), type="text",
          column.labels = c('Pooled', 'Fixed', 'Random', 'GMM'))
