FROM mcr.microsoft.com/dotnet/runtime:5.0 AS base
WORKDIR /app

FROM mcr.microsoft.com/dotnet/sdk:5.0 AS build
WORKDIR /src
COPY ["Chislaki2/Chislaki2.csproj", "Chislaki2/"]
RUN dotnet restore "Chislaki2/Chislaki2.csproj"
COPY . .
WORKDIR "/src/Chislaki2"
RUN dotnet build "Chislaki2.csproj" -c Release -o /app/build

FROM build AS publish
RUN dotnet publish "Chislaki2.csproj" -c Release -o /app/publish

FROM base AS final
WORKDIR /app
COPY --from=publish /app/publish .
ENTRYPOINT ["dotnet", "Chislaki2.dll"]
