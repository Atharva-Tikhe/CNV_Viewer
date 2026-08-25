from datetime import datetime

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from db.model import Cohort


class Service:
    def __init__(self, db: AsyncSession):
        self.db = db

    async def get_cohorts(self, id):
        result = await self.db.execute(select(Cohort).where(Cohort.id == id))
        return result.scalars().all()
