#  Pruning random points to ensure minimum distance between points
## x should be a shapefile containing the randomly generated poles

prunePole <- function(x){
  
  i=1
  
  repeat( {
    buffer <- st_buffer(x[i,], buffer_size ) #  create buffer around i-th point
    offending <- x %>%  # start with the intersection of master points... 
      st_intersects(buffer, sparse = F) # ... and the buffer, as a vector
    offending[i] <- FALSE  # i-th point, the origin, is not really offending - dont excluded
    x <- x[!offending,] 
    if ( i >= nrow(x)) {
      # the end was reached; no more points to process
      break 
    } else {
      # rinse & repeat
      i <- i + 1 
    }
  } )

  return(x)
  
}

#  Generate flight path for poles and trees within GlideRadii of a 30m pole
## x should be a shapfile containing 2019 trees and poles from ONE scenario
## within the shapefile, need to have GlideRadii (fixed for poles) and Source (North or South)
## to find out max benefits from poles
## we will only examine the flight paths involving poles (Frm or To)

pathMake <- function(x){
  
  buffer <- st_buffer(x,x$GlideRadii)
  joined <- st_join(x, buffer,
                    suffix = c(".To", ".Frm"), left = FALSE) %>%
    mutate(Type = ifelse(Type.To == "Pole", "Pole", 
                         ifelse(Type.Frm == "Pole", "Pole", "Tree"))) %>%           
    filter(Type %in% 'Pole') %>%
    filter(!ID.To == ID.Frm)
  path <- st_sfc(mapply(function(a,b){st_cast(st_union(a,b),"LINESTRING")},
                        joined$Geom.Frm, joined$Geom.To,
                        SIMPLIFY=FALSE),
                 crs=3414)
  flight <- joined %>%
            st_set_geometry(path) %>%
            dplyr::select(-c(Geom.Frm, Geom.To))
  return(flight)
  
}

# set landheight

landHeight <- function(x){
  x <- x %>%
    mutate(diffElev_Max = maxElev-Elev.To) %>%
    mutate(LandHeight = ifelse(diffElev_Max>1.5, diffElev_Max, 1.5))
  return(x)
}

# Calculate H and remove invalid flight paths

flightCheck <- function(x){

  flightValid <- x %>%
    mutate(Dist = drop_units(st_length(x))) %>%
    mutate(v = Dist/(0.8781*log(Dist) + 0.046)) %>%
    mutate(H = v + LandHeight + (Elev.To - Elev.Frm)) %>%
    #If pole, we are looking at 100% of height instead of 80%
    mutate(MaxDepart = ifelse(Type.To == "Pole",
                              v + Height_M.To + (Elev.To - Elev.Frm), #rule for pole
                              v + 0.8*Height_M.To + (Elev.To - Elev.Frm))) %>% #rule for tree
    mutate(Depart = ifelse(Type.Frm == "Pole",
                           ifelse((MaxDepart<Height_M.Frm), MaxDepart, Height_M.Frm), #rule for pole
                           ifelse((MaxDepart<0.8*Height_M.Frm), MaxDepart, 0.8*Height_M.Frm))) %>% #rule for tree
    mutate(Diff_H = Depart-H ) %>%
    mutate(Diff_M = 0.8*Height_M.To - LandHeight) %>%
    filter(! Diff_H < 0) %>%
    filter(! Diff_M < 0)
  
  #removing those that connect to the same tree i.e. serving existing connectivity
  flightValid2 <- flightValid %>%
    filter(! ID.To %in% flight22$ID_To) %>%
    filter(! ID.Frm %in% flight22$ID.Frm) %>%
    filter(! ID.Frm %in% flight22$ID_To) %>%
    filter(! ID.To %in% flight22$ID_Frm) %>%
    mutate(scale = scales::rescale(Diff_H, to = c(1.0, 1.9))) #rescaling for quality
  
  return(flightValid2)
  
}

#  Calculate valid flight path and other score metrics
## x should be valid_Flight[[s]]

#change summarise() to summaris(, .groups = 'drop') to prevent warning message from showing

scoreMetrics <- function(x){
  
  repeat({
  #unadjusted values, raw counts
  ValidSum.To_Unadjusted <- x %>%
    filter(Type.To == 'Pole') %>%
    count(ID.To, MLR_Loc.Frm) %>%
    st_drop_geometry()
  
  #going from Pole to Tree
  y_Unadjusted <- x %>%
    filter(Type.Frm == 'Pole') %>%
    count(ID.Frm, MLR_Loc.To) %>%
    st_drop_geometry() %>%
    #combine to and from
    full_join(ValidSum.To_Unadjusted, by = c("ID.Frm" = "ID.To"), keep = T) %>%
    filter(! MLR_Loc.To == MLR_Loc.Frm) %>%
    rename(ID = ID.Frm) %>%
    rename(Entry = n.y) %>%
    rename(Exit = n.x) %>%
    rename(Direction = MLR_Loc.To ) %>%
    rename(TreeLoc = MLR_Loc.Frm) %>%
    dplyr::select(- ID.To) %>%
    rowwise() %>%
    mutate(UniCont = min(Entry, Exit))
  
  #using adjusted values based on Diff_H range
  ValidSum.To <- x %>%
    filter(Type.To == 'Pole') %>%
    group_by(ID.To, MLR_Loc.Frm) %>%
    summarise(adjusted = sum(scale), .groups = 'drop') %>%
    st_drop_geometry()

  y <- x %>%
    filter(Type.Frm == 'Pole') %>%
    group_by(ID.Frm, MLR_Loc.To) %>%
    summarise(adjusted = sum(scale), .groups = 'drop') %>%
    st_drop_geometry() %>%
    #combine to and from
    full_join(ValidSum.To, by = c("ID.Frm" = "ID.To"), keep = T) %>%
    filter(! MLR_Loc.To == MLR_Loc.Frm) %>%
    rename(ID = ID.Frm) %>%
    rename(Entry = adjusted.y) %>%
    rename(Exit = adjusted.x) %>%
    rename(Direction = MLR_Loc.To ) %>%
    rename(TreeLoc = MLR_Loc.Frm) %>%
    dplyr::select(- ID.To) %>%
    rowwise() %>%
    mutate(UniCont = min(Entry, Exit)) %>%
    dplyr::left_join(y_Unadjusted %>% dplyr::select(ID, Direction, TreeLoc, UniCont),
                     by = (c("ID", "Direction", "TreeLoc"))) %>%
    rename(ScoreAdj = UniCont.x) %>%
    rename(ScoreRaw = UniCont.y) %>%
    mutate(ID.PoleLoc = paste(ID, Direction)) %>%
    mutate(ID.TreeLoc = paste(ID, Direction))

  #need to check if its pole to pole, den is the entire connection valid?
  #connected pole's Direction should be same
  
  #what are the pole IDs for pole to pole connection
  #and is part of the final set of valid poles
  test <- x[(x$Type.To == x$Type.Frm) & (x$MLR_Loc.Frm != x$MLR_Loc.To), c(1,8,11)] %>%
    st_drop_geometry() %>%
    filter(ID.To %in% y$ID | ID.Frm %in% y$ID) %>%
    mutate(Destination = NA) %>%
    mutate(Dest.FrmMiddle = NA) %>%
    mutate(Check = NA)
  
  test.sameSide <- x[(x$Type.To == x$Type.Frm) & (x$MLR_Loc.Frm == x$MLR_Loc.To), c(1,8,11)] %>%
    st_drop_geometry() %>%
    filter(ID.To %in% y$ID | ID.Frm %in% y$ID) %>%
    mutate(Destination = NA) %>%
    mutate(Dest.FrmMiddle = NA) %>%
    mutate(Check = NA)
  
  if ( nrow(test) == 0 & nrow(test.sameSide) == 0 ){ break 
    } else {
  
      for ( a in 1:nrow(test)) {
        test[a,]$Destination <- ifelse(test[a,]$ID.To %in% y$ID,
                                      as.character(y[(y$ID %in% test[a,]$ID.To) & 
                                                       (y$TreeLoc %in% test[a,]$MLR_Loc.Frm), ]$Direction),
                                      'FALSE')
        test[a,]$Dest.FrmMiddle <- ifelse(test[a,]$MLR_Loc.Frm == 'Middle', 
                                          as.character(y[(y$ID %in% test[a,]$ID.Frm), ]$Direction),
                                          "NA")
        test[a,]$Check <- ifelse(test[a,]$MLR_Loc.Frm == 'Middle',
                                 as.character(test[a,]$Destination == test[a,]$Dest.FrmMiddle),
                                 ifelse(test[a,]$Destination == FALSE, FALSE, 
                                   as.character(test[a,]$Destination != test[a,]$MLR_Loc.Frm)))
      } 
  
      test <- test %>%
      replace_na(list(Check = FALSE))
      
      if(nrow(test.sameSide) == 0) { test.sameSide <- test.sameSide
        } else {
      
        for( b in 1:nrow(test.sameSide)){
          test.sameSide[a,]$Destination <- ifelse(test.sameSide[a,]$ID.To %in% y$ID,
                                                  as.character(y[(y$ID %in% test.sameSide[a,]$ID.To) & 
                                                                   (y$TreeLoc %in% test.sameSide[a,]$MLR_Loc.Frm), ]$Direction),
                                                  'FALSE')
          test.sameSide[a,]$Dest.FrmMiddle <- ifelse(test.sameSide[a,]$MLR_Loc.Frm == 'Middle', 
                                                     as.character(y[(y$ID %in% test.sameSide[a,]$ID.Frm), ]$Direction),
                                                     "NA")
          test.sameSide[a,]$Check <- ifelse(test.sameSide[a,]$MLR_Loc.Frm == 'Middle',
                                            as.character(test.sameSide[a,]$Destination == test.sameSide[a,]$Dest.FrmMiddle),
                                            ifelse(test.sameSide[a,]$Destination == FALSE, FALSE, 
                                                   as.character(test.sameSide[a,]$Destination != test.sameSide[a,]$MLR_Loc.Frm)))
        }
        
        test.sameSide <- test.sameSide %>%
          replace_na(list(Check = FALSE))
      
      }
  
      #check ID.To's direction opposite of ID.Frm's location
      #TRUE = Not equal = opposite direction --> keep
      if(all(test$Check == TRUE) & all(test.sameSide$Check == TRUE)) { break 
        } else {
          #Check = FALSE means it's a dead end --> remove from valid flight and recalculate y
          #need to adjust to only remove the invalid pole, and not all pole
          test2 <- test %>% filter(Check == FALSE)
          x <- x %>%
            filter(! (ID.To %in% test2$ID.To & ID.Frm %in% test2$ID.Frm))
            
          if(nrow(test.sameSide) == 0) { test.sameSide2 <- test.sameSide
           } else {
            test.sameSide2 <- test.sameSide %>% filter(Check == FALSE)
           }
          
          x <- x %>%
            filter(! (ID.To %in% test.sameSide2$ID.To & ID.Frm %in% test.sameSide2$ID.Frm))
          
          #Check = TRUE but there are extra invalid flight even tho parts of it is valid
          #e.g. There are paths entering from South and North into Pole A
          #but Pole A to Pole B only allow pathway from South to North
          #which means that those coming into Pole A from North is useless
          #unless it can exit from Pole A on its own
          test3 <- test %>% rename(InvalidEntry = Dest.FrmMiddle)
          
          for ( a in 1:nrow(test)) {
            test3[a,]$Destination <- ifelse(test3[a,]$ID.To %in% y$ID,
                                            as.character(y[(y$ID %in% test3[a,]$ID.To) & 
                                                             (y$TreeLoc %in% test3[a,]$MLR_Loc.Frm), ]$Direction),
                                            'FALSE')
            test3[a,]$InvalidEntry <- ifelse(test3[a,]$MLR_Loc.Frm %in% 'South', 'North', 'South')
          }
          test3 <- test3 %>% filter(Check == TRUE)
          x <- x %>%
            filter(! (ID.To %in% test3$ID.Frm & MLR_Loc.Frm %in% test3$InvalidEntry))
        }
    }
  
})
  
  y <- y %>%
    filter(! ((ID %in% test$ID.Frm) & (TreeLoc %in% test$Destination) & (Direction %in% 'Middle')))
  
return(y)
  
}

#Creating a new scenario based on original
#shift one pole at a time and test conditions

evoShift <- function(x){
  t = 1
  Coord <- list()
  repeat({
    
    #Create a new set of Coord for scenario
    for (i in 1:nrow(x)) {
        p <- 1
        #shifting one at a time and testing if it remains within planting verge
        repeat({
          pairs <- cbind(sample(c(-step_size:step_size),1), sample(c(-step_size:step_size),1))
          Coord[i] <- x[i,]$geometry + pairs
          if (!is.na(as.integer(st_intersects(Coord[[i]],mlr_2019, sparse = TRUE)))){
            break
            } else{
              p <- p+1
          }
        })
      }

    #Update set of geometry for evo_shift(x)
    shift <- st_set_geometry(x, st_geometry(st_sfc(Coord)))%>%
      st_set_crs(st_crs(x))
    
      #Though unlikely, but still check that it fulfill minimum spacing of 25m
      if (as.numeric(min(st_distance(shift)
                          [st_distance(shift) != min(st_distance(shift))])) > 25){
        break 
        } else {
        # rinse & repeat
        t <- t + 1 
      }
  })

  return(shift)
}

## Previous evoShift function (without restraining movement to within road)
evoShift.old <- function(x){
  
  t <- 1
  
  repeat( {
    
    print(paste("try", t, sep = " "))
    
    #random sets of movement
    val1 <- sample(c(-step_size:step_size), samples, replace = TRUE)
    val2 <- sample(c(-step_size:step_size), samples, replace = TRUE)
    pairs <- cbind(val1, val2)
    
    #all poles shift acc to varying directions as specified in pairs
    Coord <- list()
    for (i in 1:nrow(x)) {
      Coord[i] <- x[i,]$geometry + pairs[i,]
    }
    
    shift <- st_set_geometry(x, st_geometry(st_sfc(Coord)))%>%
      st_set_crs(st_crs(x))
    
    if ( as.numeric(min(st_distance(shift)
                        [st_distance(shift) != min(st_distance(shift))])) > 25){
      break 
    } else {
      # rinse & repeat
      t <- t + 1 
    }
  }
  )
  
  return(shift)
  
}
