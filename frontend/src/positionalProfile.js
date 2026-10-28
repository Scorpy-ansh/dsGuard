function finiteNumber(value) {

  if (value === null || value === undefined || value === "") return null;

  const number = Number(value);

  return Number.isFinite(number) ? number : null;

}



function contiguousSegments(points) {

  const segments = [];

  let segment = [];



  points.forEach((point) => {

    if (point.value === null) {

      if (segment.length) segments.push(segment);

      segment = [];

    } else {

      segment.push(point);

    }

  });



  if (segment.length) segments.push(segment);

  return segments;

}



export function positionalProfileData(profile) {

  const targetLength = finiteNumber(profile?.target_length);

  if (!Number.isFinite(targetLength) || targetLength < 1) {

    return { available: false, reason: "Positional target length is unavailable." };

  }



  const windows = (Array.isArray(profile?.windows) ? profile.windows : [])

    .map((row) => {

      const start = finiteNumber(row.start);

      const end = finiteNumber(row.end);

      return {

        start,

        end,

        potency: finiteNumber(row.potency),

        row,

      };

    })

    .filter((row) => row.start !== null && row.end !== null && row.end >= row.start);



  const profileRows = (Array.isArray(profile?.profiles) ? profile.profiles : [])

    .map((row) => ({

      position: finiteNumber(row.position_1based),

      row,

    }))

    .filter((item) => item.position !== null && item.position >= 1);



  const species = profileRows.length

    ? Object.keys(profileRows[0].row).filter((key) => key !== "position_1based")

    : [];

  const profilePoints = (name) =>

    profileRows.map(({ position, row }) => ({

      position,

      value: finiteNumber(row[name]),

    }));



  if (!windows.length && !profileRows.length) {

    return { available: false, reason: "This run did not record valid positional data." };

  }



  return {

    available: true,

    targetLength,

    windows,

    species,

    profilePoints,

    potencyPoints: windows.map(({ start, potency }) => ({

      position: start,

      value: potency,

    })),

  };

}



export function numericRange(points) {

  const values = points.map((point) => point.value).filter((value) => value !== null);

  if (!values.length) return null;

  return { min: Math.min(...values), max: Math.max(...values) };

}



export { contiguousSegments };
